
#include <iostream>
using namespace std; 
#include <math.h>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fstream>

const double c_light = 29979245800;
const double pi = atan(1.0) * 4;
const double b = 0.942137; // this is the beta of the tau in c.m 
//vit is printed if the initial data are changed


double random_uniform(double min, double max)
{
    return min + (double)(rand() / (double)(RAND_MAX) * (max - min));
}


double distribution(double x) {
    return 1+pow(x,2);
}

double distributionTWO(double x) {
    return 1+ x*x+(1-b*b)*(1- x*x);
}

double  max_finder(double f(double), double n, double a,double b) {
    double max[2];
    for (int o = 0; o < n; o++) {
        max[1] = f(a);
        max[0] = a;
        if (f(a + o * (b - a)/n) > max[1]) {
            max[1] = f(a + o * (b - a)/n );
            max[0] = a + o * (b - a) / n;
        }
    }
    return max[1];
}

double  min_finder(double f(double), double n, double a, double b) {
    static double min[2];
    for (int l = 0; l < n; l++) {
        min[1] = f(a);
        min[0] = a;
        if (f(a + l * (b - a) / n) < min[1]) {
        min[1] = f(a + l * (b - a) / n);
        min[0] = a + l * (b - a) / n;
        }
    }
    return min[1];
}

double ACC_REG(double f(double), double n, double a, double b) { // f(double) is the distribution
    double X, Y;
    X = random_uniform(a, b);
    Y = max_finder(f, n, a, b); //find the maximun of f() in the interval [a,b] in n steps
    double r = random_uniform(0, 1);
    if ((r * Y) < f(X)) {
        return X;
    }
    else {
        while ((r * Y) > f(X)) {
            X = random_uniform(a, b);
            r = random_uniform(0, 1);
        }
        return X;
    }
}

double Ec(double mA, double mX, double mB)
{
    return (mX * mX + mA * mA - mB * mB) / (2 * mA);
}
double p_x(double p, double costheta, double phi)
{
    double sintheta = sqrt(1 - pow(costheta, 2));
    return p * sintheta * cos(phi);
}

double p_y(double p, double costheta, double phi)
{
    double sintheta = sqrt(1 - pow(costheta, 2));
    return p * sintheta * sin(phi);
}

double p_z(double p, double costheta)
{
    return p * costheta;
}

double flight_length(double p[4], double tau) {
    double v, gamma, t,l,l_mean;
    v = sqrt(pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2)) / p[0];
    gamma = 1 / sqrt(1 - pow(v, 2));
    t = gamma * tau;
    l_mean = v * t * c_light;
    double u = random_uniform(0, 1);
        l = -l_mean * log(1 - u);
        while (l > 500) {
            u = random_uniform(0, 1);
            l = -l_mean * log(1 - u);
    }
    return l;
}



int main()
{
    double ma = 1.77, mx = 1.11, mb = 0.139, E_cm = 10.56, E_el = 9.0;
    double E_pos, beta, gamma;
    int control, n,p;
    double f[3];
    srand( (unsigned) time(0));

    E_pos = pow(E_cm,2) / (4 * E_el); //positron energy
    beta = (E_el-E_pos) / (E_el+E_pos); // beta CDM
    gamma = 1 / sqrt(1 - pow(beta,2)); //gamma CDM

    cout << " E_pos ";
    cout << E_pos;
    cout << " beta";
    cout << beta;
    cout << "gamma";
    cout << gamma;

    cout << " Number of iteration: ";
    cin >> n;
    cout << " Type 0 uniform distribution, 1 to use 1+cos^2, 2 for the total one: ";
    cin >> control;
    cout << " If the initial data are changed type 1, otherwise type 0 ";
    cin >> p;

    cout << " Number of steps in which the interval where will be found max and min will be divided: ";
    cin >> f[0];
    cout << " Left extreme of the interval:";
    cin >> f[1];
    cout << " Right extreme of the interval: ";
    cin >> f[2];

    double pamu[4];
    double pa_lab[4];
    double pxmu[4];
    double px_lab[4];
    double px_com[4];
    double va[4];
    double vac[4];
    double tau_x = 2.631 * pow(10, -10);
    double tau_tau=2.903 * pow(10, -13);
    ofstream outFile;
    outFile.open("output.txt");

    for (int i = 0; i < n; i++) {
 
        pamu[0] = E_cm / 2;
        double pa = sqrt(pamu[0] * pamu[0] - ma * ma);
        if (p == 1) {
            cout << pa / pamu[0];
            return 0;
        }
        else if (p==0){
            cout << "La variabile b iniziale va bene cosi'.";
            p = 4;
        }

        double c_theta_cm;

        if (control == 0) {
            c_theta_cm = random_uniform(-1, 1);
        }
        else if(control==1){

            c_theta_cm = ACC_REG(distribution,f[0],f[1],f[2]);
        }
        else {
            c_theta_cm = ACC_REG(distributionTWO, f[0], f[1], f[2]);
        }

        outFile << c_theta_cm;
        outFile << " ";
  

        double phi_cm = random_uniform(0, 2 * pi);

        outFile << phi_cm;
        outFile << " ";

        double c_theta_x = random_uniform(-1, 1);

        outFile << c_theta_x;
        outFile << " ";

        double phi_x = random_uniform(0, 2 * pi);

        outFile << phi_x;
        outFile << " ";

        //4 - momenta di A in COM frame



        pamu[1] = p_x(pa, c_theta_cm, phi_cm);
        pamu[2] = p_y(pa, c_theta_cm, phi_cm);
        pamu[3] = p_z(pa, c_theta_cm);

        outFile << pamu[0];
        outFile << " ";
        outFile << pamu[1];
        outFile << " ";
        outFile << pamu[2];
        outFile << " ";
        outFile << pamu[3];
        outFile << " ";

            //4 - momenta di A in LAB frame

        pa_lab[0] = gamma * pamu[0] + beta * gamma * pamu[3];
        pa_lab[1] = pamu[1];
        pa_lab[2] = pamu[2];
        pa_lab[3] = (gamma * (pamu[3] + beta * pamu[0]));

        outFile << pa_lab[0];
        outFile << " ";
        outFile << pa_lab[1];
        outFile << " ";
        outFile << pa_lab[2];
        outFile << " ";
        outFile << pa_lab[3];
        outFile << " ";

            //4 - momenta di X(A->X + B) in A frame

        pxmu[0] = Ec(ma, mx, mb);
        double px = sqrt(pxmu[0] * pxmu[0] - mx * mx);
        pxmu[1] = p_x(px, c_theta_x, phi_x);
        pxmu[2] = p_y(px, c_theta_x, phi_x);
        pxmu[3] = p_z(px, c_theta_x);

        outFile << pxmu[0];
        outFile << " ";
        outFile << pxmu[1];
        outFile << " ";
        outFile << pxmu[2];
        outFile << " ";
        outFile << pxmu[3];
        outFile << " ";
 
            //boost to take X back in LAB frame
      
        va[1] = pa_lab[1] / (pa_lab[0]);
        va[2] = pa_lab[2] / (pa_lab[0]);
        va[3] = pa_lab[3] / (pa_lab[0]);

        double vq = va[1] * va[1] + va[2] * va[2] + va[3] * va[3];
        double gammat = 1 / sqrt(1 - vq);

        px_lab[0] = gammat * (pxmu[0] + va[1] * pxmu[1] + va[2] * pxmu[2] + va[3] * pxmu[3]);
        px_lab[1] = gammat * va[1] * pxmu[0] + pxmu[1] * (1 + (gammat - 1) * (va[1] * va[1]) / (vq)) + pxmu[2] * (gammat - 1) * (va[1] * va[2]) / (vq)+pxmu[3] * (gammat - 1) * (va[1] * va[3]) / (vq);
        px_lab[2] = gammat * va[2] * pxmu[0] + pxmu[2] * (1 + (gammat - 1) * (va[2] * va[2]) / (vq)) + pxmu[1] * (gammat - 1) * (va[1] * va[2]) / (vq)+pxmu[3] * (gammat - 1) * (va[3] * va[2]) / (vq);
        px_lab[3] = gammat * va[3] * pxmu[0] + pxmu[3] * (1 + (gammat - 1) * (va[3] * va[3]) / (vq)) + pxmu[2] * (gammat - 1) * (va[3] * va[2]) / (vq)+pxmu[1] * (gammat - 1) * (va[1] * va[3]) / (vq);

        outFile << px_lab[0];
        outFile << " ";
        outFile << px_lab[1];
        outFile << " ";
        outFile << px_lab[2];
        outFile << " ";
        outFile << px_lab[3];
        outFile << " ";

        //boost to take X back in COM frame

        vac[1] = pamu[1] / (pamu[0]);
        vac[2] = pamu[2] / (pamu[0]);
        vac[3] = pamu[3] / (pamu[0]);

        double vqc = vac[1] * vac[1] + vac[2] * vac[2] + vac[3] * vac[3];
        double gammac = 1 / sqrt(1 - vqc);

        px_com[0] = gammac * (pxmu[0] + vac[1] * pxmu[1] + vac[2] * pxmu[2] + vac[3] * pxmu[3]);
        px_com[1] = gammac * vac[1] * pxmu[0] + pxmu[1] * (1 + (gammac - 1) * (vac[1] * vac[1]) / (vqc)) + pxmu[2] * (gammac - 1) * (vac[1] * vac[2]) / (vqc)+pxmu[3] * (gammac - 1) * (vac[1] * vac[3]) / (vqc);
        px_com[2] = gammac * vac[2] * pxmu[0] + pxmu[2] * (1 + (gammac - 1) * (vac[2] * vac[2]) / (vqc)) + pxmu[1] * (gammac - 1) * (vac[1] * vac[2]) / (vqc)+pxmu[3] * (gammac - 1) * (vac[3] * vac[2]) / (vqc);
        px_com[3] = gammac * vac[3] * pxmu[0] + pxmu[3] * (1 + (gammac - 1) * (vac[3] * vac[3]) / (vqc)) + pxmu[2] * (gammac - 1) * (vac[3] * vac[2]) / (vqc)+pxmu[1] * (gammac - 1) * (vac[1] * vac[3]) / (vqc);


        //module of p x in cm

        double modpx;

        modpx = sqrt(px_com[1] * px_com[1] + px_com[2] * px_com[2] + px_com[3] * px_com[3]);
        outFile << modpx;
        outFile << " ";
        double fl;
//flight length
        fl = flight_length(px_lab, tau_x);
        outFile << fl;
        outFile << " ";
        fl = flight_length(pa_lab, tau_tau);
        outFile << fl;
        outFile << "\n";




    }

    outFile.close();
}

