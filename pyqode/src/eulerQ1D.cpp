
#include <Eigen/Eigen>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <filesystem>
#include <vector>

using namespace std;
using namespace Eigen;

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, 0, ", ", "\n");

VectorXd pEoS(VectorXd E, VectorXd rho){
    double gamma = 1.4;
    return E.array()*rho.array()*(gamma-1.0);
}

MatrixXd cleanMatrix(MatrixXd M){
    vector<Matrix<double,1,3> > buffer; //not sure of the type name for the rows
    VectorXd zero(M.cols()); //or appropriate comparable type
    for(int i = 0; i < M.rows(); i++){ //note: possibly a function call each time
        if ((M.row(i).array() != 0.0).any() || i == 0){
            buffer.push_back(M.row(i));
        }
    }
    MatrixXd return_value(buffer.size(), 3);
    for(int i = buffer.size(); i --> 0;){
        
        return_value.row(i) = buffer[i];
    }
    return return_value;
}

string hashTags(int n){
    /* Preenche terminal com n hashtags */
    string output;

    for (int i = 0; i < n; i++)
        output.append("#");
    
    return output;
}

string nulls(int n){
    /* Preenche terminal com n vazios */
    string output;

    for (int i = 0; i < n; i++)
        output.append(" ");
    
    return output;
}

string eigenVersion(){
    /* Returna a versão do Eigen da pasta ./eigen usada na compilação*/
    string output;
    output.append(to_string(EIGEN_WORLD_VERSION));
    output.append(".");
    output.append(to_string(EIGEN_MAJOR_VERSION));
    output.append("."); 
    output.append(to_string(EIGEN_MINOR_VERSION));
    return output;
}

void initPrint(ostream &of){
    /* Print inicial */

    of << hashTags(80) << endl;
    of << nulls(27) 
        <<"-*- Q1D Euler Solver -*-" <<  endl; 
    of << "Eigen " << eigenVersion() << endl;
    of << "Allan Moreira de Carvalho" <<  endl;
    of << hashTags(80) << endl;
}

istream& getl(istream& file, string& line){
    while (getline(file, line))
    {
        if (line[0] == '#' || 
            (line[0] == '/' && line[1] == '/') ||
            line[0] == '%' ||
            line[0] == ';' 
            ) 
            continue;
        else
            break;
    }
    return file;
}

void saveTxt(MatrixXd m, string filename){
        ofstream ofile(filename);
        ofile << m.format(CSVFormat);
        ofile.close();
}

MatrixXd readMatrix(string fileToOpen)
{

	vector<double> matrixEntries;
	ifstream matrixDataFile(fileToOpen);
	string matrixRowString;
	string matrixEntry;
	int matrixRowNumber = 0;
	while (getline(matrixDataFile, matrixRowString)) // here we read a row by row of matrixDataFile and store every line into the string variable matrixRowString
	{
		stringstream matrixRowStringStream(matrixRowString); //convert matrixRowString that is a string to a stream variable.

		while (getline(matrixRowStringStream, matrixEntry, ',')) // here we read pieces of the stream matrixRowStringStream until every comma, and store the resulting character into the matrixEntry
		{
			matrixEntries.push_back(stod(matrixEntry));   //here we convert the string to double and fill in the row vector storing all the matrix entries
		}
		matrixRowNumber++; //update the column numbers
	}
	return Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);
}

VectorXd cellCenters(VectorXd xn){
    int N = xn.size();
    double dxi = xn(1) - xn(0);
    double dxf = xn(N-1) - xn(N-2);
    VectorXd xc(N+1);

    for (int i = 1; i < xc.size()-1; i++)
    {
        xc(i) = (xn(i)+xn(i-1))/2.0;
    }
    xc(0) = xn(0) - (xn(1)-xn(0))/2.0;
    xc(N) = xn(N-1) + (xn(N-1)-xn(N-2))/2.0;

    return xc;
}

class fluid{
    public:
        double R, gamma;

        void setR(double R){
            this->R=R;
        }
        void setGamma(double gamma){
            this->gamma=gamma;
        }
        VectorXd rhoEoS(VectorXd p, double T, double R){
            return p.array()/(R*T);
        }
        VectorXd cEoS(VectorXd p, VectorXd rho){
            return sqrt(gamma*p.array()/rho.array());
        }
        VectorXd eEoS(VectorXd p, VectorXd rho){
            return p.array()/(gamma-1.0)/rho.array();
        }
        VectorXd pEoS(VectorXd e, VectorXd rho){
            return e.array()*rho.array()*(gamma-1.0);;
        }
};



class domain{
    public:
        VectorXd xn, xc, Sn, Sc, dx, dS;
        
        void setX(VectorXd xn){
            this->xn = xn;
            xc = cellCenters(xn);
            dx = diff(xn);
            dx(0) = dx(1);
            dx(dx.size()-1) = dx(dx.size()-2); 
        }
        VectorXd diff(VectorXd x){
            VectorXd d(x.size()+1);
            for (int i = 1; i < d.size()-1; i++)
            {
                d(i) = x(i)-x(i-1);
            }
            return d;
        }
        void setS(VectorXd Sn){
            this->Sn = Sn;
            Sc = cellCenters(Sn);
            dS = diff(Sn);
            dS(0) = dS(1);
            dS(dS.size()-1) = dS(dS.size()-2);
        }
};

class euler{
    public:
        domain d;
        fluid f;
        string dim, tscheme, fscheme, dttype, terminalfile, outpath;

        int itmax, itprint;
        double p0in, T0in, pb, Min, Tin, pin,
               CFL, tol, R, gamma,
               TRef, pRef, RRef, cRef, rhoRef, LRef;

        VectorXd rho, u, e, p, c, M, T, dt, E;

        MatrixXd U, Unew, F, Q, res, resit, convF, dU, Fmais, Fmenos;

        ofstream outfile;

        void setup(domain d, double p0in, double T0in, double Min, double pb, fluid f){
            this->d=d;
            this->f=f;
            this->p0in=p0in;
            this->T0in=T0in;
            this->pb=pb;
            this->Min=Min;
            this->R = f.R;
            this->gamma = f.gamma;
        }

        void setupFile(string filename){

            ifstream file(filename);
            string line;
            string xnFile;
            string SnFile;

            // setup domain
            getl(file, line);
            xnFile = line.c_str();

            getl(file, line);
            SnFile = line.c_str();
            
            VectorXd xn = readMatrix(xnFile).col(0);
            VectorXd Sn = readMatrix(SnFile).col(0);

            d.setX(xn);
            d.setS(Sn);

            // setup boundary conditions
            getl(file, line);
            p0in = stod(line);

            getl(file, line);
            T0in = stod(line);

            getl(file, line);
            Min = stod(line);

            getl(file, line);
            pb = stod(line);

            // setup fluid
            getl(file, line);
            R = stod(line);

            getl(file, line);
            gamma = stod(line);

            f.setR(R);
            f.setGamma(gamma);

            // setup solver
            getl(file, line);
            itmax = stoi(line);

            getl(file, line);
            itprint = stoi(line);

            getl(file, line);
            CFL = stod(line);

            getl(file, line);
            tol = stod(line);

            getl(file, line);
            tscheme = line.c_str();

            getl(file, line);
            fscheme =line.c_str();

            getl(file, line);
            dttype = line.c_str();

            getl(file, line);
            dim = line.c_str();

            getl(file, line);
            outpath = line.c_str();
            terminalfile = outpath + "terminal";
            
            filesystem::create_directory(outpath);
            outfile.open(terminalfile);

            file.close();
            // print setup file

            initPrint(cout);
            initPrint(outfile);
            cout << nulls(27) 
               <<"-*- Setup Information -*-" <<  endl;
            cout << hashTags(80) << endl;
            outfile << nulls(27) 
               <<"-*- Setup Information -*-" <<  endl;
            outfile << hashTags(80) << endl;

            printSetup(filename);
        }

        void printSetup(string filename){
            ifstream file(filename);
            if (file.is_open()) {
                string line;
                while (getline(file, line)) {
                    cout << line.c_str() << endl;
                    outfile << line.c_str() << endl;

                }
                file.close();
            }
        }

        //void solveOld(int itmax, int itprint, double CFL, double tol, 
        //           string tscheme, string fscheme, string dttype, string dim){
        void solve(){
            
            cout << hashTags(80) << endl;
            cout << nulls(27) 
                 <<"-*- Start Solving -*-" <<  endl;  
            cout << hashTags(80) << endl;

            outfile << hashTags(80) << endl;
            outfile << nulls(27) 
                    <<"-*- Start Solving -*-" <<  endl;
            outfile << hashTags(80) << endl;

            setIC(dim);

            resit = MatrixXd::Zero(itmax+1,3);
            res = MatrixXd::Zero(U.rows(),U.cols());
            dU = MatrixXd::Zero(U.rows(),U.cols());
            Unew = MatrixXd::Zero(U.rows(),U.cols());
            this->CFL = CFL;
            this->itmax = itmax;
            this->tol = tol;
            double maxres = 2*tol;

            
            //setDt(CFL, dttype);
            //updatePrimitives();
            //setFlux(fscheme);
            //setSource();

            cout <<  "\t \t \t     Maximum Residue" << endl;
            outfile <<  "\t \t \t     Maximum Residue" << endl;
                    
            cout <<  "Iteration \tMass\t\tMomentum\tEnergy" << endl;
            outfile <<  "Iteration \tMass\t\tMomentum\tEnergy" << endl;

            for (int it = 1; it < itmax+1 && maxres > tol ; it++)
            {
                updatePrimitives();
                setDt(CFL, dttype);
                setSource();
                setFlux(fscheme);

                Unew = U;
                if (tscheme == "Euler"){
                    for(int i=1; i < U.rows()-1; i++){
                        res.row(i) =  Q.row(i) - (F.row(i)*d.Sn(i) - F.row(i-1)*d.Sn(i-1))/(d.dx(i)*d.Sc(i));
                        Unew.row(i) = Unew.row(i) + res.row(i)*dt(i);
                    }
                }

                if (tscheme == "RK4"){
                    for (int rki = 1; rki < 5; rki++)
                    {
                        for(int i=1; i < U.rows()-1; i++){
                            //res.row(i) = (F.row(i)*d.Sn(i) - F.row(i-1)*d.Sn(i-1))/d.dx(i)/d.Sc(i) - Q.row(i);
                            res.row(i) = Q.row(i) - (F.row(i)*d.Sn(i) - F.row(i-1)*d.Sn(i-1))/(d.dx(i)*d.Sc(i));
                            Unew.row(i) = Unew.row(i) + res.row(i)*dt(i)/(5.0 - double(rki));
                        }
                    }
                }
        
                dU = Unew-U;
                U = Unew;

                updatePrimitives();
                setBC();  
                
                U = Unew + dU;

                
                resit.row(it-1) = dU.block(1,0,U.rows()-2,U.cols()).array().abs().colwise().maxCoeff();
                

                if (it%itprint == 0){
                    cout << it << "\t\t" << resit(it-1,0) << "\t" << resit(it-1,1) << "\t" << resit(it-1,2) << endl;
                    outfile << it << "\t\t" << resit(it-1,0) << "\t" << resit(it-1,1) << "\t" << resit(it-1,2) << endl;
                }

                maxres = resit.row(it-1).maxCoeff();
            }

            resit = cleanMatrix(resit);

            cout << hashTags(80) << endl;
            cout << nulls(27) 
               <<"-*- Saving Results -*-" <<  endl; 
            cout << hashTags(80) << endl;


            outfile << hashTags(80) << endl;
            outfile << nulls(27) 
               <<"-*- Saving Results -*-" <<  endl; 
            outfile << hashTags(80) << endl;
            
        }

        void reconstruct(){
            for(int i=1; i < U.rows()-1; i++){
                U.row(i) = U.row(i-1) + (d.xc(i))*( U.row(i+1) - U.row(i-1) )/ (d.xc(i+1)- d.xc(i-1));
                //U.row(i) = U.row(i);    
            }
        }

        VectorXd linearInterpolation(VectorXd y){
            for(int i=1; i < y.size()-1; i++){
                //y(i) = y(i-1) + (d.xc(i))*( y(i+1) - y(i-1) )/ (d.xc(i+1)- d.xc(i-1));
                y(i) = (y(i+1)+y(i-1))/2.0;
            }
            return y;
        }

        void reconstructPrimitives(){
            rho = linearInterpolation(rho);
            u = linearInterpolation(u);
            e = linearInterpolation(e);
            p = linearInterpolation(p);
            c = linearInterpolation(c);
            M = linearInterpolation(M);
            T = linearInterpolation(T);
        }
        
        void updatePrimitives(){
            rho = U.col(0).array();
            u = U.col(1).array()/U.col(0).array();
            e = U.col(2).array();
            p = (U.col(2).array()-(U.col(1).array()*U.col(1).array()/U.col(0).array())/2.0)*(gamma-1.0);
            T = p.array()/rho.array()/R;
            c = sqrt(gamma*p.array()/rho.array());
            M = u.array()/c.array();

            //E = U.col(2).array()/U.col(0).array() - 0.5*((U.col(1)*U.col(1))/(U.col(0)*U.col(0)))
            //p = EoS(E,rho);
            //T = EoS(E,rho);
            //c = EoS(E,rho);
            //M = EoS(E,rho);
        }

        void setBC(){
            inletBC();
            outletBC();
        }

        void setBC2(){
            inletBC2();
            outletBC2();
        }
        
        
        void inletBC2(){
            // Algo de errado nao esta certo
            Min = (M(0)+M(1))/2.0;
            if (Min < 1){
                double uav = (u(0)+u(1))/2.0;
                double cav = (c(0)+c(1))/2.0;
                double rhoav = (rho(0)+rho(1))/2.0;
                
                double dudx = (u(1) - u(0))/d.dx(1);
                double dpdx = (p(1) - p(0))/d.dx(1);

                double R3 = -(uav - cav)*(dudx - dpdx/(rhoav*cav)) + uav*cav*d.dS(1)/d.dx(1)/d.Sn(1);

                double dpdu = -p0in*gamma*Min/cav*pow(1.0+(gamma-1.0)*pow(Min,2.0)/2.0,(1.0-2.0*gamma)/(gamma-1.0));

                double du = R3*dt(1)/(1.0-(dpdu/(rhoav*cav)));

                double unew = u(0) + du;

                double a2 = 2.0*gamma*T0in/(gamma+1.0);

                double Tnew  = T0in*(1.0-((gamma-1.0)/(gamma+1.0))*(pow(unew,2.0)/(a2*a2)));
                double pnew  = p0in*pow(Tnew/T0in,gamma/(gamma-1.0));
                double rhonew = pnew/Tnew/R;
                double enew  = pnew/(gamma-1.0)+ rhonew*pow(unew,2.0)/2.0;
                double cnew    = sqrt(R*gamma*Tnew);
                double Mnew = unew/cnew;
    
                rho(0)  = rhonew;
                u(0)    = unew;
                c(0)    = cnew;
                p(0)    = pnew;
                T(0)    = Tnew;
                e(0)    = enew;
                M(0)    = Mnew;

                
                dU(0,0) = rhonew - U(0,0); 
                dU(0,1) = rhonew*unew - U(0,1);
                dU(0,2) = enew - U(0,2);
            }
        }
        
        
        void inletBC(){
            
            double rhog = U(0,0);
            double ug = U(0,1)/rhog;
            double pg = (gamma-1.0)*(U(0,2)-(U(0,1)*U(0,1)/U(0,0))/2.0);
            double cg = sqrt(gamma*pg/rhog);

            double rhod = U(1,0);
            double ud = U(1,1)/rhod;
            double pd = (gamma-1.0)*(U(1,2)-(U(1,1)*U(1,1)/U(1,0))/2.0);
            double cd = sqrt(gamma*pd/rhod);
            
            double Min = ug/cg;
            if (Min < 1.0){
                double gm1 = gamma-1;
                double gp1 = gamma+1;
                //double a2 = 2.0*gamma*R*T(0)/(gamma+1.0);
                double a2 = 2.0*gamma*R*T0in/(gamma+1.0);
                double dpdu = p0in*(gamma/gm1)*pow((1.0-(gm1/gp1)*ug*ug/a2),(1.0/gm1))*(-2.0*(gm1/gp1)*ug/a2);
                double dtdx = dt(1)/d.dx(1);
                double eig = ((ug-cg + ud-cd)/2.0)*dtdx;
                double dpdx = pd-pg;
                double dudx = ud-ug;

                double du = -eig*(dpdx-rhog*cg*dudx)/(dpdu-rhog*cg);

                double ugnew = ug+du;
                double Tgnew = T0in*(1.0-(gm1/gp1)*ugnew*ugnew/a2);
                double pgnew = p0in*pow((Tgnew/T0in),(gamma/gm1));
                double rhognew = pgnew/(R*Tgnew);
                double egnew = pgnew/gm1 + rhognew*0.5*ug*ug;
                
                //double cv = R/(gamma-1.0);
                //double egnew = rhognew*(cv*Tgnew+0.5*ugnew*ugnew);

                dU(0,0) = rhognew - U(0,0); 
                dU(0,1) = rhognew*ugnew - U(0,1);
                dU(0,2) = egnew - U(0,2);

            } 
        }
        
        void outletBC2(){
            // Algo de errado nao esta certo
            int N = U.rows();

            double uav = (u(N-2) + u(N-1))/2.0;
            double cav = (c(N-2) + c(N-1))/2.0;
            double rhoav = (rho(N-2) + rho(N-1))/2.0;
            
            double drhodx = (rho(N-1)-rho(N-2))/d.dx(N-2);
            double dudx = (u(N-1) - u(N-2))/d.dx(N-2);
            double dpdx = (p(N-1) - p(N-2))/d.dx(N-2);
                
            double R1 = uav*dpdx/pow(cav,2.0)-uav*drhodx;
            double R2 = -(uav + cav)*(dudx + dpdx/(rhoav*cav)) - uav*cav*d.dS(N-2)/d.dx(N-2)/d.Sn(N-2);
            double R3 = -(uav - cav)*(dudx - dpdx/(rhoav*cav)) + uav*cav*d.dS(N-2)/d.dx(N-2)/d.Sn(N-2);

            double Mout = (M(N-2)+M(N-1))/2.0;
            if (Mout <= 1){          
                double drho = dt(N-2)*R1;
                double du = dt(N-2)*R2;
            
                double pnew = p(N-1);    
                double rhonew  = rho(N-1) + drho;
                double unew    = u(N-1) + du;
                double enew    = pnew/(gamma-1.0) + rhonew*pow(unew,2.0)/2.0;
                double Tnew    = pnew/rhonew/R;
                double cnew    = sqrt(R*gamma*Tnew);
                double Mnew = unew/cnew;
    
                rho(N-1)  = rhonew;
                u(N-1)    = unew;
                c(N-1)    = cnew;
                p(N-1)    = pnew;
                T(N-1)    = Tnew;
                e(N-1)    = enew;
                M(N-1)    = Mnew;

                dU(N-1,0) = rhonew - U(N-1,0);
                dU(N-1,1) = rhonew*unew - U(N-1,1);
                dU(N-1,2) = enew - U(N-1,2);
            }
            else if (Mout > 1){
                double du   = (R2 + R3)*dt(N-2)/2.0;
                double dp   = rhoav*cav*(R2 - R3)*dt(N-2)/2.0;
                double drho = R1*dt(N-2) + dp/pow(cav,2.0); 
                
                double pnew = p(N-1) + dp;
                double rhonew  = rho(N-1) + drho;
                double unew    = u(N-1) + du;
                double enew    = pnew/(gamma-1.0) + rhonew*pow(unew,2.0)/2.0;   
                double Tnew    = pnew/rhonew/R;
                double cnew    = sqrt(R*gamma*Tnew);
                double Mnew = unew/cnew;
    
                rho(N-1)  = rhonew;
                u(N-1)    = unew;
                c(N-1)    = cnew;
                p(N-1)    = pnew;
                T(N-1)    = Tnew;
                e(N-1)    = enew;
                M(N-1)    = Mnew;

                dU(N-1,0) = rhonew - U(N-1,0);
                dU(N-1,1) = rhonew*unew - U(N-1,1);
                dU(N-1,2) = enew - U(N-1,2);
            } 
        }
        
        
        
        void outletBC(){
            double gm1 = gamma-1;
            double gp1 = gamma+1;
            int N = U.rows();
            
            double rhog = U(N-1,0);
            double ug = U(N-1,1)/rhog;
            double pg = (gamma-1.0)*(U(N-1,2)-(U(N-1,1)*U(N-1,1)/U(N-1,0))/2.0);
            double cg = sqrt(gamma*pg/rhog);

            double rhod = U(N-2,0);
            double ud = U(N-2,1)/rhod;
            double pd = (gamma-1.0)*(U(N-2,2)-(U(N-2,1)*U(N-2,1)/U(N-2,0))/2.0);
            double cd = sqrt(gamma*pd/rhod);

            double dtdx = dt(N-2)/d.dx(N-2);
            double avu = 0.5*(ud+ug);
            double avc = 0.5*(cd+cg);

            double eig0 = avu*dtdx;
            double eig1 = (avu+avc)*dtdx;
            double eig2 = (avu-avc)*dtdx;
            
            double dpdx = pg-pd;
            double dudx = ug-ud;

            double R0 = -eig0*((rhog-rhod)-dpdx/(cg*cg));
            double R1 = -eig1*(dpdx+rhog*cg*dudx);
            double R2 = -eig2*(dpdx-rhog*cg*dudx);

            double Mout = avu/avc;;
            double dp = 0;
            if (Mout > 1.0)
                dp = 0.5*(R1+R2);

            double drho = R0+dp/(cg*cg);
            double du = (R1-dp)/(rhog*cg);

            double ugnew = ug + du;
            double rhognew = rhog + drho;
            double pgnew = pg + dp;
            double Tgnew = pgnew/(rhognew*R);
            double egnew = pgnew/gm1 + rhognew*0.5*ug*ug;
            
            //double cv = R/(gamma-1.0);
            //double egnew = rhognew*(cv*Tgnew+0.5*ugnew*ugnew);

            dU(N-1,0) = rhognew - U(N-1,0);
            dU(N-1,1) = rhognew*ugnew - U(N-1,1);
            dU(N-1,2) = egnew - U(N-1,2);
        }
        
        void setFlux(string scheme){
            if (scheme == "Roe")
                setFluxRoe();
            else if (scheme == "Basic")
                setFluxDefault();
            else if (scheme == "DefaultV")
                setFluxDefaultVector();
            else if (scheme == "Split")
                setFluxSplit();
            else if (scheme == "SplitTVD")
                setFluxSplitTVD();
            else if (scheme == "VanLeer")
                setFluxVanLeer();
            else if (scheme == "AUSM") // works flawless
                setFluxAUSM();
            else if (scheme == "LaxFriedrichs")
                setFluxLaxFriedrichs(); 
            else if (scheme == "StegerWarming")
                setFluxStegerWarming();
        }

        void setFluxStegerWarming(){
            VectorXd Fp = VectorXd::Zero(3);
            VectorXd Fm = VectorXd::Zero(3);

            F = MatrixXd::Zero(U.rows(), U.cols());

            for (int i = 0; i < F.rows()-1; i++)
            {
                double ML = M(i);
                double maxML = max(0.0, ML);
                double cL = c(i);
                double rhoL = rho(i);

                double MR = M(i+1);
                double minMR = min(0.0, MR);
                double cR = c(i+1);
                double rhoR = rho(i+1);

                
                if(ML*ML > 1.0){
                    Fp(0)=rhoL*cL*maxML;
                    Fp(1)=Fp(0)*cL/(gamma*ML)*(gamma*ML*ML+1.0)*maxML;
                    Fp(2)=Fp(0)*cL*cL*0.5*(2.0/(gamma-1.0)+ML*ML)*maxML;
                }else{
                    double X = 0.5*rhoL*cL/gamma*(ML+1.0);
                    double Y = rhoL*cL*(gamma-1.0)/gamma;

                    Fp(0)=X+Y*maxML;
                    Fp(1)=X*cL*(ML+1.0)+Y*cL*ML*maxML;
                    Fp(2)=X*0.5*cL*cL*((ML+1.0)*(ML+1.0)+(3.0-gamma)/(gamma-1))+Y*0.5*cL*cL*ML*ML*maxML;
                }

                if(MR*MR > 1.0){
                    Fm(0)=rhoR*cR*minMR;
                    Fm(1)=Fm(0)*cR/(gamma*MR)*(gamma*MR*MR+1.0)*minMR;
                    Fm(2)=Fm(0)*cR*cR*0.5*(2.0/(gamma-1.0)+MR*MR)*minMR;
                }else{
                    double X = 0.5*rhoR*cR/gamma*(MR-1.0);
                    double Y = rhoR*cR*(gamma-1.0)/gamma;

                    Fm(0)=X+Y*minMR;
                    Fm(1)=X*cR*(MR-1.0)+Y*cR*MR*minMR;
                    Fm(2)=X*0.5*cR*cR*((MR-1.0)*(MR-1.0)+(3.0-gamma)/(gamma-1))+Y*0.5*cR*cR*MR*MR*minMR;
                }
                F.row(i) = Fp + Fm;
            }
        }

        void setFluxLaxFriedrichs(){
            convFlux();

            F = MatrixXd::Zero(U.rows(), U.cols());

            for (int i = 0; i < F.rows()-1; i++)
            {
                F.row(i) = 0.5*(d.dx(i)/dt(i)*(U.row(i)-U.row(i+1)) + convF.row(i) + convF.row(i+1) );
            }
        }
        
        void setFluxSplit(){
         
            double eps=1e-2;
            MatrixXd Pinv(3,3), P(3,3);
            MatrixXd Lmais = MatrixXd::Zero(3,3);
            MatrixXd Lmenos = MatrixXd::Zero(3,3);
            MatrixXd Amais = MatrixXd::Zero(3,3);
            MatrixXd Amenos = MatrixXd::Zero(3,3);
            MatrixXd Fmais = MatrixXd::Zero(U.rows(),U.cols());
            MatrixXd Fmenos = MatrixXd::Zero(U.rows(),U.cols());

            //vector<Matrix<double,3,3> > Amais(U.rows());
            //vector<Matrix<double,3,3> > Amenos(U.rows());
            

            F = MatrixXd::Zero(U.rows(), U.cols());

            for (int i = 0; i < F.rows(); i++)
            {
                // split flux method Steger Warming
                double eps = 1e-5;

                Pinv(0,0) = 1.0-(gamma-1.0)/2.0*pow(u(i),2.0)/(pow(c(i),2.0));
                Pinv(0,1) = (gamma-1.0)*u(i)/(pow(c(i),2.0));
                Pinv(0,2) = -(gamma-1.0)/(pow(c(i),2.0));
                Pinv(1,0) = ((gamma-1.0)*pow(u(i),2.0)/2.0-u(i)*c(i))/(rho(i)*c(i));
                Pinv(1,1) = (c(i)-(gamma-1.0)*u(i))/(rho(i)*c(i));
                Pinv(1,2) = (gamma-1.0)/(rho(i)*c(i));
                Pinv(2,0) = -((gamma-1.0)*pow(u(i),2)/2+u(i)*c(i))/(rho(i)*c(i));
                Pinv(2,1) = (c(i)+(gamma-1.0)*u(i))/(rho(i)*c(i));
                Pinv(2,2) = -(gamma-1.0)/(rho(i)*c(i));

                P(0,0) = 1.0;
                P(0,1) = rho(i)/(2.0*c(i));
                P(0,2) = -rho(i)/(2.0*c(i));
                P(1,0) = u(i);
                P(1,1) = rho(i)*(u(i)+c(i))/(2.0*c(i));
                P(1,2) = -rho(i)*(u(i)-c(i))/(2.0*c(i));
                P(2,0) = (pow(u(i),2.0))/2.0;
                P(2,1) = (rho(i)/(2.0*c(i)))*((pow(u(i),2.0)/2.0) + u(i)*c(i) + pow(c(i),2.0)/(gamma-1.0));
                P(2,2) = -(rho(i)/(2.0*c(i)))*((pow(u(i),2.0)/2.0) - u(i)*c(i) + pow(c(i),2.0)/(gamma-1.0));


                double l1 = u(i);
                double l2 = u(i)+c(i);
                double l3 = u(i)-c(i);

                Lmais << (l1 + sqrt(pow(l1,2.0)+pow(eps,2.0)))/2.0, 0.0,                0.0,
                        0.0,          (l2 + sqrt(pow(l2,2.0)+pow(eps,2.0)))/2.0,       0.0,
                        0.0,          0.0,       (l3 + sqrt(pow(l3,2.0)+pow(eps,2.0)))/2.0;

                

                Lmenos << (l1 - sqrt(pow(l1,2.0)+pow(eps,2.0)))/2.0, 0.0,                0.0,
                        0.0,          (l2 - sqrt(pow(l2,2.0)+pow(eps,2.0)))/2.0,       0.0,
                        0.0,          0.0,       (l3 - sqrt(pow(l3,2.0)+pow(eps,2.0)))/2.0;
                
                /*double lambda1 = u(i);
                double lambda2 = u(i) + c(i);
                double lambda3 = u(i) - c(i);
                
                Lmais(0,0) = (lambda1 + sqrt(pow(lambda1,2)+pow(eps,2)))/2;
                Lmais(1,1) = (lambda2 + sqrt(pow(lambda2,2)+pow(eps,2)))/2;
                Lmais(2,2) = (lambda3 + sqrt(pow(lambda3,2)+pow(eps,2)))/2;
                    
                Lmenos(0,0) = (lambda1 - sqrt(pow(lambda1,2)+pow(eps,2)))/2;
                Lmenos(1,1) = (lambda2 - sqrt(pow(lambda2,2)+pow(eps,2)))/2;
                Lmenos(2,2) = (lambda3 - sqrt(pow(lambda3,2)+pow(eps,2)))/2;*/
                
                Amais  = P*Lmais*Pinv;      
                Amenos = P*Lmenos*Pinv;
                
                Fmais.row(i)  = Amais*U.row(i).transpose();
                Fmenos.row(i) = Amenos*U.row(i).transpose();    

            }
            
            for (int i = 0; i < F.rows()-1; i++)
            {
                F.row(i) = Fmais.row(i) + Fmenos.row(i+1);
            }
          
        }

        void setFluxAUSM(){
            VectorXd Fmais = VectorXd::Zero(3);
            VectorXd Fmenos = VectorXd::Zero(3);

            F = MatrixXd::Zero(U.rows(), U.cols());

            for (int i = 0; i < F.rows()-1; i++)
            {
                double ML = M(i);
                double PL = p(i);
                double cL = c(i);
                double rhoL = rho(i);
                double uL = u(i);
                double HL = (e(i) + p(i))/rho(i);

                double MR = M(i+1);
                double PR = p(i+1);
                double cR = c(i+1);
                double rhoR = rho(i+1);
                double uR = u(i+1);
                double HR= (e(i+1) + p(i+1))/rho(i+1);

                double Mp = 0.0;
                double Pp = 0.0;
                if(ML <= -1){
                    //nothing
                }
                if ( ML < 1.0){
                    Mp = 0.25*(ML+1.0)*(ML+1.0);
                    Pp = 0.25*PL*(ML+1.0)*(ML+1.0)*(2.0-ML);
                }
                else{
                    Mp = ML;
                    Pp = PL;
                }

                double Mm = 0;
                double Pm = 0;
                if (MR <=-1){
                    Mm=MR;
                    Pm=PR;
                }
                else if ( MR < 1.0){
                    Mm = -0.25*(MR-1.0)*(MR-1.0);
                    Pm = 0.25*PR*(1.0-MR)*(1.0-MR)*(2.0+MR);
                }

                Fmais(0)=max(0.0,Mp+Mm)*rhoL*cL;
                Fmais(1)=max(0.0,Mp+Mm)*rhoL*cL*uL + Pp;
                Fmais(2)=max(0.0,Mp+Mm)*rhoL*cL*HL;
                
                Fmenos(0)=min(0.0,Mp+Mm)*rhoR*cR;
                Fmenos(1)=min(0.0,Mp+Mm)*rhoR*cR*uR + Pm;
                Fmenos(2)=min(0.0,Mp+Mm)*rhoR*cR*HR;

                F.row(i) = Fmais + Fmenos;
            }
        }

        void convFlux(){
            convF = MatrixXd::Zero(U.rows(), U.cols());

            convF << rho, 
                     rho.array()*u.array()*u.array()+p.array(),
                     (e.array() + p.array())*u.array();
        }

        void setFluxSplitTVD(){
         
            double eps=1e-2;
            MatrixXd Pinv(3,3), P(3,3);
            MatrixXd Lmais = MatrixXd::Zero(3,3);
            MatrixXd Lmenos = MatrixXd::Zero(3,3);
            MatrixXd Amais = MatrixXd::Zero(3,3);
            MatrixXd Amenos = MatrixXd::Zero(3,3);
            MatrixXd Fmais = MatrixXd::Zero(U.rows(),U.cols());
            MatrixXd Fmenos = MatrixXd::Zero(U.rows(),U.cols());

            //vector<Matrix<double,3,3> > Amais(U.rows());
            //vector<Matrix<double,3,3> > Amenos(U.rows());
            

            F = MatrixXd::Zero(U.rows(), U.cols());

            for (int i = 0; i < F.rows(); i++)
            {
                // split flux method Steger Warming
                double eps = 1e-5;

                Pinv(0,0) = 1.0-(gamma-1.0)/2.0*pow(u(i),2.0)/(pow(c(i),2.0));
                Pinv(0,1) = (gamma-1.0)*u(i)/(pow(c(i),2.0));
                Pinv(0,2) = -(gamma-1.0)/(pow(c(i),2.0));
                Pinv(1,0) = ((gamma-1.0)*pow(u(i),2.0)/2.0-u(i)*c(i))/(rho(i)*c(i));
                Pinv(1,1) = (c(i)-(gamma-1.0)*u(i))/(rho(i)*c(i));
                Pinv(1,2) = (gamma-1.0)/(rho(i)*c(i));
                Pinv(2,0) = -((gamma-1.0)*pow(u(i),2)/2+u(i)*c(i))/(rho(i)*c(i));
                Pinv(2,1) = (c(i)+(gamma-1.0)*u(i))/(rho(i)*c(i));
                Pinv(2,2) = -(gamma-1.0)/(rho(i)*c(i));

                P(0,0) = 1.0;
                P(0,1) = rho(i)/(2.0*c(i));
                P(0,2) = -rho(i)/(2.0*c(i));
                P(1,0) = u(i);
                P(1,1) = rho(i)*(u(i)+c(i))/(2.0*c(i));
                P(1,2) = -rho(i)*(u(i)-c(i))/(2.0*c(i));
                P(2,0) = (pow(u(i),2.0))/2.0;
                P(2,1) = (rho(i)/(2.0*c(i)))*((pow(u(i),2.0)/2.0) + u(i)*c(i) + pow(c(i),2.0)/(gamma-1.0));
                P(2,2) = -(rho(i)/(2.0*c(i)))*((pow(u(i),2.0)/2.0) - u(i)*c(i) + pow(c(i),2.0)/(gamma-1.0));


                double l1 = u(i);
                double l2 = u(i)+c(i);
                double l3 = u(i)-c(i);
                eps = l3*0.3;

                Lmais << (l1 + sqrt(pow(l1,2.0)+pow(eps,2.0)))/2.0, 0.0,                0.0,
                        0.0,          (l2 + sqrt(pow(l2,2.0)+pow(eps,2.0)))/2.0,       0.0,
                        0.0,          0.0,       (l3 + sqrt(pow(l3,2.0)+pow(eps,2.0)))/2.0;

                

                Lmenos << (l1 - sqrt(pow(l1,2.0)+pow(eps,2.0)))/2.0, 0.0,                0.0,
                        0.0,          (l2 - sqrt(pow(l2,2.0)+pow(eps,2.0)))/2.0,       0.0,
                        0.0,          0.0,       (l3 - sqrt(pow(l3,2.0)+pow(eps,2.0)))/2.0;
                
                /*double lambda1 = u(i);
                double lambda2 = u(i) + c(i);
                double lambda3 = u(i) - c(i);
                
                Lmais(0,0) = (lambda1 + sqrt(pow(lambda1,2)+pow(eps,2)))/2;
                Lmais(1,1) = (lambda2 + sqrt(pow(lambda2,2)+pow(eps,2)))/2;
                Lmais(2,2) = (lambda3 + sqrt(pow(lambda3,2)+pow(eps,2)))/2;
                    
                Lmenos(0,0) = (lambda1 - sqrt(pow(lambda1,2)+pow(eps,2)))/2;
                Lmenos(1,1) = (lambda2 - sqrt(pow(lambda2,2)+pow(eps,2)))/2;
                Lmenos(2,2) = (lambda3 - sqrt(pow(lambda3,2)+pow(eps,2)))/2;*/
                
                Amais  = P*Lmais*Pinv;      
                Amenos = P*Lmenos*Pinv;
                
                Fmais.row(i)  = Amais*U.row(i).transpose();
                Fmenos.row(i) = Amenos*U.row(i).transpose();    

            }
            
            for (int i = 0; i < F.rows()-1; i++)
            {
                F.row(i) = Fmais.row(i) + Fmenos.row(i+1);
            }
          
        }

        void setFluxDefault(){
            double eps = 0.2;
            convF = MatrixXd::Zero(U.rows(), U.cols());
            F = MatrixXd::Zero(U.rows(), U.cols());

            convF << rho, 
                     rho.array()*u.array()*u.array()+p.array(),
                     (e.array() + p.array())*u.array();

                              
            for (int i = 0; i < convF.rows()-1 ; i++)
            {
                double uav = 0.5*(u(i)+u(i+1));
                double cav = 0.5*(c(i)+c(i+1));
                double lambda = uav + cav;

                F.row(i) = 0.5*(convF.row(i+1)+convF.row(i)) - eps*lambda*(U.row(i+1)-U.row(i));
            }      
            
        }

        void setFluxDefaultVector(){

            double eps = 1e-5;
            convF = MatrixXd::Zero(U.rows(), U.cols());
            F = MatrixXd::Zero(U.rows(), U.cols());

            convF << rho, 
                     rho.array()*u.array()*u.array()+p.array(),
                     (e.array() + p.array())*u.array();

            for (int i = 0; i < F.rows()-1; i++)
            {
                double uav = 0.5*(u(i)+u(i+1));
                double cav = 0.5*(c(i)+c(i+1));

                double lambda1 = uav;
                double lambda2 = uav + cav;
                double lambda3 = uav - cav;
                
                F(i,0) = 0.5*(convF(i+1,0)+convF(i,0)) - eps*lambda1*(U(i+1,0)-U(i,0));
                F(i,1) = 0.5*(convF(i+1,1)+convF(i,1)) - eps*lambda2*(U(i+1,1)-U(i,1));
                F(i,2) = 0.5*(convF(i+1,2)+convF(i,2)) - eps*lambda3*(U(i+1,2)-U(i,2));
                
            }
        }

        void setFluxRoe(){
            F = MatrixXd::Zero(U.rows(), U.cols());
            VectorXd dV= VectorXd::Zero(3); 
            VectorXd lambda = VectorXd::Zero(3);
            MatrixXd P = MatrixXd::Zero(3, 3);

            convFlux();
            
            
            for (int i = 0; i < F.rows()-1; i++)
            {
                double HL = (e(i) + p(i))/rho(i);
                double HR = (e(i+1) + p(i+1))/rho(i+1);

                double R = sqrt(rho(i+1)/rho(i));
                double ry = rho(i)*R;
                double uy = (R*u(i+1)+u(i))/(R+1.0);
                double hy = (R*HR+HL)/(R+1.0);
                double ay = sqrt((gamma-1.0)*(hy-0.5*uy*uy));

                double drho = rho(i+1)-rho(i);
                double du = u(i+1)-u(i);
                double dp = p(i+1)-p(i);

                dV << 0.5*(dp-ry*ay*du)/ay/ay,
                      -(dp/ay/ay-drho),
                      0.5*(dp+ry*ay*du)/ay/ay;
                
                lambda << abs(uy-ay), 
                          abs(uy),
                          abs(uy+ay);

                
                double Da = max(0.0, 4.0*((u(i+1)-c(i+1))-(u(i)-c(i))));
                if ( lambda(0) < 0.5*Da )
                    lambda(0) = lambda(0)*lambda(0)/Da + 0.25*Da;
                Da = max(0.0, 4.0*((u(i+1)+c(i+1))-(u(i)+c(i))));
                if ( lambda(2) < 0.5*Da )
                    lambda(2) = lambda(2)*lambda(2)/Da + 0.25*Da;

                P << 1.0,      1.0,     1.0 ,
                     uy-ay,    uy,        uy+ay   ,
                     hy-ay*uy, 0.5*uy*uy, hy+ay*uy;

                F.row(i) = 0.5*(convF.row(i)+convF.row(i+1));

                
                for (int j = 0; j < F.cols(); j++)
                {
                    for (int k = 0; k < F.cols(); k++)
                    {
                        F(i,j) = F(i,j) -0.5*lambda(k)*dV(k)*P(j,k);
                    }
                }
            }
        }

        void setFluxOldRoe(){
            VectorXd h0 = VectorXd::Zero(rho.size());
            h0 =  gamma/(gamma-1.0)*p.array()/rho.array()+0.5*pow(u.array(),2.0);

            MatrixXd Froe = MatrixXd::Zero(U.rows(), U.cols());
            convF = MatrixXd::Zero(U.rows(), U.cols());
            F = MatrixXd::Zero(U.rows(), U.cols());

            convF << rho, 
                     rho.array()*u.array()*u.array()+p.array(),
                     (e.array() + p.array())*u.array();
            
            VectorXd dU = VectorXd::Zero(3); 

            MatrixXd Pinv = MatrixXd::Zero(3, 3);
            MatrixXd P = MatrixXd::Zero(3, 3);
            MatrixXd lambda = MatrixXd::Zero(3, 3);
            MatrixXd A = MatrixXd::Zero(3, 3);

            for (int i = 0; i < Froe.rows()-1; i++)
            {
                double R = sqrt(rho(i+1)/rho(i));
                double ry = rho(i)*R;
                double uy = (R*u(i+1)+u(i))/(R+1.0);
                double hy = (R*h0(i+1)+h0(i))/(R+1.0);
                double ay = sqrt((gamma-1.0)*(hy-0.5*uy*uy));

                double alp1 = (gamma-1.0)*uy*uy/(2.0*ay*ay);
                double alp2 = (gamma-1.0)/(ay*ay);

                dU = U.row(i+1)-U.row(i);

                Pinv << 0.5*(alp1+uy/ay), -0.5*(alp2*uy+1.0/ay), alp2/2.0,
                        1.0-alp1          , alp2*uy            , -alp2 ,
                        0.5*(alp1-uy/ay), -0.5*(alp2*uy-1.0/ay), alp2/2.0;

                P << 1.0,        1.0,         1.0 ,
                     uy-ay,    uy,        uy+ay   ,
                     hy-ay*uy, 0.5*uy*uy, hy+ay*uy;

                lambda << abs(uy-ay),   0.0,       0.0,
                          0.0,          abs(uy), 0.0,
                          0.0,          0.0,       abs(uy+ay);

                A = P*lambda*Pinv;
                
                Froe.row(i) = A*dU;
                
                F.row(i) = 0.5*(convF.row(i)+convF.row(i+1))-0.5*Froe.row(i);
            }
        }

        void setFluxRoeOld(){
            VectorXd Q = VectorXd::Zero(3);
            VectorXd h0 = VectorXd::Zero(rho.size());
            h0 =  gamma/(gamma-1.0)*p.array()/rho.array()+0.5*pow(u.array(),2.0);

            MatrixXd R = MatrixXd::Zero(3,3);
            MatrixXd L = MatrixXd::Zero(3,3);
            MatrixXd A = MatrixXd::Zero(3,3);

            VectorXd dU = VectorXd::Zero(3); 

            MatrixXd convF = MatrixXd::Zero(U.rows(), U.cols());
            convF = MatrixXd::Zero(U.rows(), U.cols());

            F = MatrixXd::Zero(U.rows(), U.cols());

            MatrixXd Rinv = MatrixXd::Zero(3,3);

            convF << rho, 
                     rho.array()*u.array()*u.array()+p.array(),
                     (e.array() + p.array())*u.array();

            for(int i = 0; i < F.rows()-1; i++){
                Q(0) = 0.5*(sqrt(rho(i)) +sqrt(rho(i+1)));
                Q(1) = 0.5*(sqrt(rho(i))*u(i) + sqrt(rho(i+1))*u(i+1));
                Q(2) = 0.5*(sqrt(rho(i))*h0(i) + sqrt(rho(i+1))*h0(i+1));

                double uy = Q(1)/Q(0);
                double hy = Q(2)/Q(0);
                double ay = sqrt((gamma-1)*(hy-0.5*uy*uy));

                R << 1,         1,          1,
                     uy - ay,   uy,         uy + ay,
                     hy - uy*ay, 0.5*uy*uy, hy+uy*ay;

                L << abs(uy-ay), 0,       0,
                     0,          abs(uy), 0,
                     0,          0,       abs(uy+ay);

                Rinv = R.inverse();
                A = R*L*Rinv;
                dU = U.row(i+1)-U.row(i);
                
                F.row(i) = 0.5*(convF.row(i) + convF.row(i+1)) - 0.5*(A*dU).transpose();
            }
            
        }

        void setFluxVanLeer(){
            F = MatrixXd::Zero(U.rows(), U.cols());
            
            VectorXd Fp = VectorXd::Zero(3);
            VectorXd Fm = VectorXd::Zero(3);

            for (int i = 0; i < F.rows()-1; i++)
            {
                double ML = M(i);
                double cL = c(i);
                double rhoL = rho(i);

                double MR = M(i+1);
                double cR = c(i+1);
                double rhoR = rho(i+1);

                Fp(0)=0.25*rhoL*cL*(ML+1.0)*(ML+1.0);
                Fp(1)=Fp(0)*2.0*cL*(1.0+0.5*(gamma-1.0)*ML)/gamma;
                Fp(2)=Fp(0)*2.0*cL*cL*pow(1.0+0.5*(gamma-1.0)*ML, 2.0)/(gamma*gamma-1.0);

                Fm(0)=-0.25*rhoR*cR*(MR-1.0)*(MR-1.0);
                Fm(1)=Fm(0)*2.0*cR*(-1.0+0.5*(gamma-1.0)*MR)/gamma;
                Fm(2)=Fm(0)*2.0*cR*cR*pow(1.0-0.5*(gamma-1.0)*MR, 2.0)/(gamma*gamma-1.0);

                F.row(i) = Fp + Fm;
            }
            
        }

        void setFluxVanLeerOld(){
            F = MatrixXd::Zero(U.rows(), U.cols());
            Fmais = MatrixXd::Zero(U.rows(), U.cols());
            Fmenos = MatrixXd::Zero(U.rows(), U.cols());

            for (int i = 0; i < F.rows()-1; i++)
            {
                // Toro1997
                Fmais.row(i) << 1.0, 
                                2.0*c(i)/gamma*((gamma-1.0)/2.0*M(i)+1.0),
                                2.0*c(i)*c(i)/(gamma*gamma-1.0)*pow((gamma-1.0)/2.0*M(i)+1.0, 2.0);
                Fmais.row(i) = 1.0/4.0*rho(i)*c(i)*pow(1.0 + M(i),2) * Fmais.row(i);

                Fmenos.row(i) << 1.0, 
                                 2.0*c(i)/gamma*((gamma-1.0)/2.0*M(i)-1.0),
                                 2.0*c(i)*c(i)/(gamma*gamma-1.0)*pow((gamma-1.0)/2.0*M(i)-1.0, 2.0);
                Fmenos.row(i) = -1.0/4.0*rho(i)*c(i)*pow(1.0-M(i),2)*Fmenos.row(i);             
            }

            for (int i = 0; i < F.rows()-1; i++)
            {
                if (M(i) >= 1)
                    F.row(i) = Fmenos.row(i);
                else
                    F.row(i) = Fmenos.row(i) + Fmenos.row(i+1);
            }

        }
        
        void setDtGlobal(double CFL){
            dt = VectorXd::Ones(d.xc.size());
            int N = u.size();
            double dt_temp;
            double dt_ = CFL*d.dx(1)/abs(u(1)+c(1));
            for(int i=1; i < N-1; i++){
                if(u(0)!=0 && c(0)!=0){
                dt_temp = CFL*d.dx(i)/abs(u(i)+c(i));
                    if ( dt_temp < dt_){
                        dt_ = dt_temp;
                    }
                }
            }
            dt = dt*dt_;
        }
        
        void setDt(double CFL, string type){
            if (type == "Local")
                setDtLocal(CFL);
            else if (type == "Global")
                setDtGlobal(CFL);
            
        }

        void setDtLocal(double CFL){
            dt = VectorXd::Zero(d.xc.size());
            dt = (CFL*d.dx.array())/abs(u.array()+c.array());
        }

        void setIC(string type){
            dim = type;

            Tin = T0in/(1.0+(gamma-1.0)/2.0*Min*Min);
            pin = p0in*pow((1.0+(gamma-1.0)/2.0*Min*Min),(-gamma/(gamma-1.0)));
            
            if (dim == "Dimensionless"){
                // // adimensionalize the Domain too;
                LRef = d.xn(d.xn.size()-1) - d.xn(0);
                d.xn = d.xn/LRef;
                d.xc = d.xc/LRef;
                d.dx = d.dx/LRef;

                TRef = T0in;
                pRef = p0in;
                RRef = R;    
                R = 1.0;
                rhoRef = pRef/RRef/TRef;
                cRef = sqrt(gamma*pRef/rhoRef);

                pin = pin/(rhoRef*cRef*cRef);
                Tin = Tin/TRef;
                pb = pb/rhoRef/cRef/cRef;

                T0in = 1.0;
                p0in = p0in/(rhoRef*cRef*cRef);
                
            }

            

            
            p = VectorXd::LinSpaced(d.xc.size(), pin, pb);
            rho = f.rhoEoS(p,Tin,R);
            c = f.cEoS(p,rho);
            u = c.array()*Min;
            M = u.array()/c.array();
            //double cv = R/(gamma-1.0);
            //e = rho.array()*(cv*Tin+0.5*u.array()*u.array());
            e = p.array()/(gamma-1.0)+ rho.array() * pow(u.array(),2.0)/2.0;
            T = p.array()/rho.array()/R;

            U = MatrixXd::Zero(d.xc.size(), 3);
            U << rho, rho.array()*u.array(), e.array();

            Q = MatrixXd::Zero(U.rows(),U.cols());
            setSource();

        }
        void setSource(){
            Q.col(1) = p.array()/d.Sc.array() * d.dS.array()/d.dx.array();
            //Q.col(1) = p.array()* d.dS.array();
        }

        void save(string path){
            filesystem::create_directory(path);
            // Dimensionalization
            if (dim == "Dimensionless"){
                p = p.array()*rhoRef*cRef*cRef;
                T = T.array()*TRef;
                rho = p.array()/RRef/T.array();
                e = e.array()*rhoRef*cRef*cRef;
                c = sqrt(gamma*p.array()/rho.array());
                u = M.array()*c.array();
                d.xc = d.xc.array()*LRef;
            }
        
            saveTxt(p                 ,path+"p.txt");
            saveTxt(M                 ,path+"M.txt");
            saveTxt(T                 ,path+"T.txt");
            saveTxt(rho               ,path+"rho.txt");
            saveTxt(e                 ,path+"e.txt");
            saveTxt(c                 ,path+"c.txt");
            saveTxt(u                 ,path+"u.txt");
            saveTxt(d.xc              ,path+"x.txt");
            saveTxt(d.Sc              ,path+"S.txt");
        
            
            saveTxt(resit.col(0),path+"resrho.txt");
            saveTxt(resit.col(1),path+"resrhou.txt");
            saveTxt(resit.col(2),path+"rese.txt");
        }
};


VectorXd S(VectorXd x, double a, double b, double c){
    //return VectorXd::Ones(x.size());
    return 1.0 - a*pow(sin(M_PI*pow(x.array(),b)),c);
}


VectorXd SYeom(VectorXd x, double xt){
    VectorXd s = VectorXd::Zero(x.size());

    for (int i = 0; i < x.size(); i++)
    {
        if ( x(i) < xt)
            s(i) = 2.5 + 3.0*(x(i)/xt - 1.5)*(x(i)/xt)*(x(i)/xt);
        if ( x(i) >= xt)
            s(i) = 3.5 - (x(i)/xt)*(6.0 - 4.5*(x(i)/xt) + (x(i)/xt)*(x(i)/xt));
    }

    return s;
    
}

int main(int argc, char const *argv[])
{
    
    euler e;
    
    e.setupFile(argv[1]);
    e.solve();
    e.save(e.outpath);

    return 0;
}
