#include<iostream>
#include<vector>
#include<math.h>
#include<cmath>
#include<climits>
#include<fstream>
#include<iomanip>
using namespace std;
double largest(int a, int b, int c)
{
    if(a>b) {
        if(a>=c) {
            return a;
        } else{
            return c;
        }
    } else{
        if(b >= c) {
            return b;
        } else {
            return c;
        }
    }
}
double smallest(int a, int b, int c)
{
    if (a <= b && a <= c)
        return a;
 
    else if (b <= a && b <= c)
        return b;
 
    else
        return c;
}
double gibbs(double zl,double zg,double A,double B)
{
    double del1=2.414;
    double del2=-0.414;
    double difference=((zg-zl)+log((zl-B)/(zg-B))-((A/(B*(del2-del1)))*log(((zl+(del1*B))/(zl+(del2*B)))*((zg+(del2*B))/(zg+(del1*B))))));
    if (difference>=0)
    {
        return zl;
    }
    else{
        return zg;
    }
}
double ze(double A,double B)
{
     double a, b, c, d;
    a=1;
    b=B-1;
    c=(A-(2*B)-(3*B*B));
    d=-((A*B)-(B*B)-(B*B*B));
    b= b/a;
    c= c/a;
    d= d/a;
    
    double discriminant, qi, ri, du, s, t, te, r;
    qi = (3.0*c - (b*b))/9.0;
    ri = -(27.0*d) + b*(9.0*c - 2.0*(b*b));
    ri =ri/ 54.0;
    discriminant = qi*qi*qi + ri*ri;
    te = (b/3.0);
    
    double x1, x2, x3;
    double x2_imag, x3_imag;
    string x2_imag_s, x3_imag_s;
    if (discriminant > 0)
    {
        s = ri + sqrt(discriminant);
        s = s<0 ? -cbrt(-s) : cbrt(s);
        t = ri - sqrt(discriminant);
        t = t<0 ? -cbrt(-t) : cbrt(t);
        x1 = -te + s + t;
        te += (s + t)/2.0;
        x3 = x2 = -te;
        te= sqrt(3.0)*(-t + s)/2;
        x2_imag = te;
        x3_imag = -te;
        x2_imag_s =  " + "+ to_string(x2_imag) + "i";
        x3_imag_s =  " - "+ to_string(x2_imag) + "i";
        return x1;
    } 
    else if (discriminant == 0)  
    { 
        x3_imag = x2_imag = 0;
        r = ri<0 ? -cbrt(-ri) : cbrt(ri);
        x1 = -te + 2.0*r;
        x3 = x2 = -(r + te);
        double zl=smallest(x1,x2,x3);
        double zg=largest(x1,x2,x3);
        return gibbs(zl,zg,A,B);
    }
    else
    {
        x3_imag = x2_imag = 0;
        qi = -qi;
        du = qi*qi*qi;
        du = acos(ri/sqrt(du));
        r = 2.0*sqrt(qi);
        x1 = -te + r*cos(du/3.0);
        x2 = -te + r*cos((du + 2.0*M_PI)/3.0);
        x3 = -te + r*cos((du + 4.0*M_PI)/3.0);
        double zl=smallest(x1,x2,x3);
        double zg=largest(x1,x2,x3);
        return gibbs(zl,zg,A,B);
    }
}
double phie(double Z,double A,double B)
{
    double res=(Z-1)-log(Z-B)+((A/(B*2*1.414))*log((Z-(0.414*B))/(Z+(2.414*B))));
    return exp(res);
}
double hi(double delB,double temperature,double pressure,double fh,double density)
{
    double res=(1.114535*log(fh))-(0.114535*log((83.144*temperature*density)/18))+(2*density*delB);
    return exp(res);
}
double density(double temp,double pressure)
{
    temp=temp-273.15;
    double a1= 3.2891-(0.002391*temp)+(0.00028446*temp*temp)-(0.00000282*pow(temp,3))+(8.477*pow(10,-9)*pow(temp,4));
    double a2= 6.245*pow(10,-5)-(3.913*pow(10,-6)*temp)-(3.499*pow(10,-8)*temp*temp)+(7.942*pow(10,-10)*pow(temp,3))-(3.299*pow(10,-12)*pow(temp,4));
    double B= 19654.32+(147.037*temp)-(2.2155*temp*temp)+(0.010478*pow(temp,3))-(2.2789*pow(10,-5)*pow(temp,4));
    double v=(1+(0.0181597*temp))/(0.9998+(18.2249*pow(10,-3)*temp)-(7.9222*pow(10,-6)*temp*temp)-(55.4485*pow(10,-9)*pow(temp,3))+(149.7562*pow(10,-12)*pow(temp,4))-(393.2952*pow(10,-15)*pow(temp,5)));
    double res=v-((v*pressure)/(B+(a1*pressure)+(a2*pressure*pressure)));
    return 1/res;
}
double fh(double T,double Tch20,double Pch20,double den,double pressure)
{

    double Psfug = Pch20 * exp((Tch20 / T) * (-7.8595178 * (1 - T / Tch20) + 1.8440825 * pow((1 - T / Tch20), 1.5) -11.786649 * pow((1 - T / Tch20), 3) + 22.680741 * pow((1 - T / Tch20), 3.5) -15.9618719 * pow((1 - T / Tch20), 4) + 1.8012250 * pow((1 - T / Tch20), 7.5)));
    double res=Psfug*exp((18.0152*(pressure-Psfug))/(den*83.144*T));
    return res;
}
int main()
{
    double temperature,pressure,brine;
      ofstream op("output.txt");
   op<<setw(17)<<"TEMPERATURE"<<setw(17)<<"PRESSURE"<<setw(17)<<"Z      "<<setw(17)<<"FUGACITY COEFFICIENT       "<<setw(17)<<"HENRY COEFFICIENT   KCO2        "<<"     KH20"<<setw(17)<<"YH20"<<setw(17)<<"XCO2"<<setw(17)<<endl;
    ifstream gp("Input.txt");
    while ( gp >> temperature>>pressure>>brine)
    {
    double w=0.224;
    double tc=304.2;
    double tch2o=647.1;
    double pc=73.83;
    double pch2o=220.55;
    double r=0.083144;
    double k=0.37464+(1.54226*w)-(0.26992*w*w);
    double tr=temperature/tc;
    double pr=pressure/pc;
    double a=(0.45724*r*r*tc*tc)*(1+k*(1-(sqrt(tr))))*(1+k*(1-(sqrt(tr))))/pc;
    double b=0.07780*r*tc/pc;
    double A=(a*pressure)/(r*r*temperature*temperature);
    double B=(b*pressure)/(r*temperature);
    double Z=ze(A,B);
    double ph=phie(Z,A,B);
    double delB=(6.187967*sqrt(1000/temperature))-5.279063;
    double dens=density(temperature,pressure);
    double fnoth2o=fh(temperature,tch2o,pch2o,dens,pressure);
    double hco2=hi(delB,temperature,pressure,fnoth2o,dens);
    double lamda=-0.0652869+(1.6790636E-04*temperature)+(40.838951/temperature)-(3.9266518E-2*pressure/temperature)+((2.1157167E-02*pressure)/(630-temperature))+(6.5486487E-06*temperature*log(pressure));
    double epsilon=-1.144624E-2+(2.8274958E-5*temperature)+((1.3980876E-2*pressure)/temperature)-((1.4349005E-2*pressure)/(630-temperature));
    double gamma=exp((2*brine*lamda)+(2*brine*brine*epsilon));
    double kco2=(hco2*gamma)/(pressure*ph);
    double knot=pow(10,(-2.209+(3.097*pow(10,-2)*(temperature-273.15))-(1.098*pow(10,-4)*(temperature-273.15)*(temperature-273.15))+(2.048*pow(10,-7)*pow(temperature-273.15,3))));
    double kh20=(knot/(fnoth2o*pressure))*exp(((pressure-1)*18.18)/(83.14*temperature));
    double yh20=(1-(1/kco2))/((1/kh20)-(1/kco2));
    double ynormal=1/(1+yh20);
    double xi=((ynormal/kco2));
        op<<setw(17)<<temperature<<setw(17)<<pressure<<setw(17)<<Z<<setw(17)<<ph<<setw(17)<<hco2<<setw(17)<<kco2<<setw(17)<<kh20<<setw(17)<<yh20<<setw(17)<<xi<<endl; 
    }

    gp.close();
    op.close();
}