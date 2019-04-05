#include <iostream>
#include <math.h>
using namespace std;


int main(){

    const double pi = 3.14159265358979;
	double m=0.0; //month
	double Day=0.0; //Date
	double h=0.0;
	double min=0.0;
	double s=0.0;
    double y=0.0; //year
	double d=7019.92;
    cout<<"Please enter the year:  ";
	cin>>y;
	cout<<"Please enter the month:  ";
	cin>>m;
	cout<<"Please enter the day:  ";
	cin>>Day;
	cout<<"Please enter the hour (in 24 hours format):  ";
	cin>>h;
	cout<<"Please enter the minute:  ";
	cin>>min;
	cout<<"Please enter the second:  ";
	cin>>s;
	double UT;
    UT = h + min/60 + s/3600; //UT in hours + decimals
    d = 367*y - 7 * ( y + (m+9)/12) /4 + 275*m/9 + Day - 730530 ; //day
	d = d + UT/24.0;
    double ecl = 23.4393 - 3.563E-7 * d;  //obliquity of the ecliptic, i.e. the "tilt" of the Earth's axis of rotation
	//cout<< d<<endl;
    ////////////////////////////////Sun////////////////////////////////
    double N_S=0.0;  //longitude of the ascending node
    double i_S=0.0;  //inclination to the ecliptic (plane of the Earth's orbit)
    double w_S=282.9404 + 4.70935E-5 * d;  //argument of perihelion
    double a_S=1.000000;  //semi-major axis, or mean distance from Sun (AU)
    double e_S=0.016709 - 1.151E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    double M_S=356.0470 + 0.9856002585 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
	while(M_S>(360))
		M_S = M_S - 360;
	while(M_S<0)
		M_S = M_S + 360;
	double L_S= M_S + w_S; //Mean Longitude of the Sun  (Ns=0)
    ////////////////////////////////Moon////////////////////////////////
    double N_M=125.1228 - 0.0529538083 * d;  //longitude of the ascending node
	while(N_M>360)
		N_M = N_M - 360;
	while(N_M<0)
		N_M = N_M + 360;
    double i_M=5.1454;  //inclination to the ecliptic (plane of the Earth's orbit)
    double w_M=318.0634 + 0.1643573223 * d;  //argument of perihelion
	while(w_M>(360))
		w_M = w_M - 360;
	while(w_M<0)
		w_M = w_M + 360;
    double a_M=60.2666;  //semi-major axis, or mean distance from Sun (earth radii)
    double e_M=0.054900;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    double M_M=115.3654 + 13.0649929509 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_M>(360))
		M_M = M_M - 360;
	while(M_M<0)
		M_M = M_M + 360;
	double L_M = M_M + w_M + N_M ; //Mean longitude of the Moon
	double D= L_M - L_S; //Mean elongation of the Moon
	double F= L_M - N_M; //Argument of latitude for the Moon
    ////////////////////////////////Mercury////////////////////////////////
    double N_ME=48.3313 + 3.24587E-5 * d;  //longitude of the ascending node
    double i_ME=7.0047 + 5.00E-8 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    double w_ME=29.1241 + 1.01444E-5 * d;  //argument of perihelion
    double a_ME=0.387098;  //semi-major axis, or mean distance from Sun
    double e_ME=0.205635 + 5.59E-10 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    double M_ME=168.6562 + 4.0923344368 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_ME>(360))
		M_ME = M_ME - 360;
	while(M_ME<0)
		M_ME = M_ME + 360;
    ////////////////////////////////Venus////////////////////////////////
    double N_V=76.6799 + 2.46590E-5 * d;  //longitude of the ascending node
    double i_V=3.3946 + 2.75E-8 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    double w_V=54.8910 + 1.38374E-5 * d;  //argument of perihelion
    double a_V=0.723330;  //semi-major axis, or mean distance from Sun
    double e_V=0.006773 - 1.302E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    double M_V=48.0052 + 1.6021302244 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_V>(360))
		M_V = M_V - 360;
	while(M_V<0)
		M_V = M_V + 360;
    ////////////////////////////////Mars////////////////////////////////
    double N_MA=49.5574 + 2.11081E-5 * d;  //longitude of the ascending node
    double i_MA=1.8497 - 1.78E-8 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    double w_MA=286.5016 + 2.92961E-5 * d;  //argument of perihelion
    double a_MA=1.523688;  //semi-major axis, or mean distance from Sun
    double e_MA=0.093405 + 2.516E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    double M_MA=18.6021 + 0.5240207766 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_MA>(360))
		M_MA = M_MA - 360;
	while(M_MA<0)
		M_MA = M_MA + 360;
    ////////////////////////////////Jupiter////////////////////////////////
    double N_J=100.4542 + 2.76854E-5 * d;  //longitude of the ascending node
    double i_J=1.3030 - 1.557E-7 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    double w_J=273.8777 + 1.64505E-5 * d;  //argument of perihelion
    double a_J=5.20256  ;  //semi-major axis, or mean distance from Sun
    double e_J=0.048498 + 4.469E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    double M_J=19.8950 + 0.0830853001 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_J>(360))
		M_J = M_J - 360;
	while(M_J<0)
		M_J = M_J + 360;
    ////////////////////////////////Saturn////////////////////////////////
    double N_SA=113.6634 + 2.38980E-5 * d;  //longitude of the ascending node
    double i_SA=2.4886 - 1.081E-7 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    double w_SA=339.3939 + 2.97661E-5 * d;  //argument of perihelion
    double a_SA=9.55475  ;  //semi-major axis, or mean distance from Sun
    double e_SA=0.055546 - 9.499E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    double M_SA=316.9670 + 0.0334442282 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_S>(360))
		M_S = M_S - 360;
	while(M_S<0)
		M_S = M_S + 360;
    ////////////////////////////////Uranus////////////////////////////////
    double N_U=74.0005 + 1.3978E-5 * d;  //longitude of the ascending node
    double i_U=0.7733 + 1.9E-8 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    double w_U=96.6612 + 3.0565E-5 * d;  //argument of perihelion
    double a_U=19.18171 - 1.55E-8 * d;  //semi-major axis, or mean distance from Sun
    double e_U=0.047318 + 7.45E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    double M_U=142.5905 + 0.011725806 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_U>(360))
		M_U = M_U - 360;
	while(M_U<0)
		M_U = M_U + 360;
    ////////////////////////////////Neptune////////////////////////////////
    double N_N=131.7806 + 3.0173E-5 * d;  //longitude of the ascending node
    double i_N=1.7700 - 2.55E-7 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    double w_N=272.8461 - 6.027E-6 * d;  //argument of perihelion
    double a_N=30.05826 + 3.313E-8 * d;  //semi-major axis, or mean distance from Sun
    double e_N=0.008606 + 2.15E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    double M_N=260.2471 + 0.005995147 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_N>(360))
		M_N = M_N - 360;
	while(M_N<0)
		M_N = M_N + 360;
    ////////////////////////////////////////////////////////////////


    ////////////////////////////////Sun's Position////////////////////////////////
    double E_S = M_S + e_S*(180/pi) * sin(M_S*pi/180) * ( 1.0 + e_S * cos(M_S*pi/180) );  //Sun's eccentric anomaly E_S
    double xv_S = a_S * (cos(E_S*pi/180) - e_S); //NOT IMPORTANT
    double yv_S = a_S * sqrt(1.0 - e_S*e_S) * sin(E_S*pi/180); //NOT IMPORTANT
    double vs = atan2( yv_S, xv_S ); //Sun's true anomaly
    double rs = sqrt( xv_S*xv_S + yv_S*yv_S ); //Sun's distance
    double lonsun = vs*180/pi + w_S; //sun's true longitude
	while(lonsun>360)
		lonsun = lonsun - 360;
    double xs = rs * cos(lonsun*pi/180); //geocentric coordinates
    double ys = rs * sin(lonsun*pi/180); //geocentric coordinates zs=0
    double xes = xs; //equatorial coordinates
    double yes = ys * cos(ecl*pi/180); //equatorial coordinates
    double zes = ys * sin(ecl*pi/180); //equatorial coordinates
    double RA_S  = atan2( yes, xes )*180/pi; //Sun's Right Ascension
    double Dec_S = atan2( zes, sqrt(xes*xes+yes*yes) )*180/pi; //Sun's Declination
    ////////////////////////////////////////////////////////////////

    ////////////////////////////////Sidreal Time////////////////////////////////
	double GMST0 = L_S + 180; //Greenwich Mean Sidereal Time at 0h UT 
	double GMST = GMST0 + UT * 15.0; //The Greenwich Mean Sideral Time
	double LST  = GMST + 51.5; //replace 51.5 (the longitude of Tehran) by your local longitude
    ////////////////////////////////////////////////////////////////

    ////////////////////////////////Moon's Position////////////////////////////////
    double E_M = M_M + e_M*(180/pi) * sin(M_M*pi/180) * ( 1.0 + e_M * cos(M_M*pi/180) );  //Moon's eccentric anomaly E_S
    if(e_M>0.05){
        double E0_M=E_M;
        double E1_M=E_M + 2;
        while(fabs(E0_M-E1_M)>0.001){
            E0_M=E1_M;
            E1_M = E_M - ( E_M - e_M*(180/pi) * sin(E_M*pi/180) - M_M ) / ( 1 - e_M * cos(E_M*pi/180) );
        }
		E_M=E1_M;
    }
    double xv_M = a_M * (cos(E_M*pi/180) - e_M); //NOT IMPORTANT
    double yv_M = a_M * sqrt(1.0 - e_M*e_M) * sin(E_M*pi/180); //NOT IMPORTANT
    double v_M = atan2( yv_M, xv_M ); //Moon's true anomaly
    double r_M = sqrt( xv_M*xv_M + yv_M*yv_M ); //Moon's distance
    double xh_M = r_M * ( cos(N_M*pi/180) * cos(v_M+w_M*pi/180) - sin(N_M*pi/180) * sin(v_M+w_M*pi/180) * cos(i_M*pi/180) ); //Moon's heliocentric (Sun-centered) position in the ecliptic coordinate system
    double yh_M = r_M * ( sin(N_M*pi/180) * cos(v_M+w_M*pi/180) + cos(N_M*pi/180) * sin(v_M+w_M*pi/180) * cos(i_M*pi/180) ); //Moon's heliocentric (Sun-centered) position in the ecliptic coordinate system
    double zh_M = r_M * ( sin(v_M+w_M*pi/180) * sin(i_M*pi/180) ); //Moon's heliocentric (Sun-centered) position in the ecliptic coordinate system
	double lonecl_M = atan2( yh_M, xh_M )*180/pi; //Moon's ecliptic longitude
    double latecl_M = atan2( zh_M, sqrt(xh_M*xh_M+yh_M*yh_M) )*180/pi; //Moon's ecliptic latitude
	////////////////////////////////Mercury's Position////////////////////////////////
    double E_ME = M_ME + e_ME*(180/pi) * sin(M_ME*pi/180) * ( 1.0 + e_ME * cos(M_ME*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_ME>0.05){
        double E0_ME=E_ME;
        double E1_ME=E_ME;
        while(fabs(E0_ME-E1_ME)>0.001){
            E0_ME=E1_ME;
            E1_ME = E_ME - ( E_ME - e_ME*(180/pi) * sin(E_ME*pi/180) - M_ME ) / ( 1 - e_ME * cos(E_ME*pi/180) );
        }
    }
    double xv_ME = a_ME * (cos(E_ME*pi/180) - e_ME); //NOT IMPORTANT
    double yv_ME = a_ME * sqrt(1.0 - e_ME*e_ME) * sin(E_ME*pi/180); //NOT IMPORTANT
    double v_ME = atan2( yv_ME, xv_ME ); //Mercury's true anomaly
    double r_ME = sqrt( xv_ME*xv_ME + yv_ME*yv_ME ); //Mercury's distance
    double xh_ME = r_ME * ( cos(N_ME*pi/180) * cos(v_ME+w_ME*pi/180) - sin(N_ME*pi/180) * sin(v_ME+w_ME*pi/180) * cos(i_ME*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
    double yh_ME = r_ME * ( sin(N_ME*pi/180) * cos(v_ME+w_ME*pi/180) + cos(N_ME*pi/180) * sin(v_ME+w_ME*pi/180) * cos(i_ME*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
    double zh_ME = r_ME * ( sin(v_ME+w_ME*pi/180) * sin(i_ME*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	double lonecl_ME = atan2( yh_ME, xh_ME )*180/pi; //Mercury's ecliptic longitude
    double latecl_ME = atan2( zh_ME, sqrt(xh_ME*xh_ME+yh_ME*yh_ME) )*180/pi; //Mercury's ecliptic latitude
	double xg_ME = r_ME * cos(lonecl_ME*pi/180) * cos(latecl_ME*pi/180) +xs; //geocentric position
	double yg_ME = r_ME * sin(lonecl_ME*pi/180) * cos(latecl_ME*pi/180) +ys; //geocentric position
	double zg_ME = r_ME * sin(latecl_ME*pi/180); //geocentric position
	double xe_ME = xg_ME; //equatorial coordinates
    double ye_ME = yg_ME * cos(ecl*pi/180) - zg_ME * sin(ecl*pi/180); //equatorial coordinates
    double ze_ME = yg_ME * sin(ecl*pi/180) + zg_ME * cos(ecl*pi/180); //equatorial coordinates
    double RA_ME  = atan2( ye_ME, xe_ME )*180/pi; //Right Ascension
    double Dec_ME = atan2( ze_ME, sqrt(xe_ME*xe_ME+ye_ME*ye_ME) )*180/pi; //Declination
	double par_ME = (8.794/3600) / r_ME;
	double HA_ME = LST - RA_ME; //hour angle of the Moon
	double x_ME = cos(HA_ME*pi/180) * cos(Dec_ME*pi/180);
    double y_ME = sin(HA_ME*pi/180) * cos(Dec_ME*pi/180);
    double z_ME = sin(Dec_ME*pi/180);
	double xhor_ME = x_ME * sin(35*pi/180) - z_ME * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
    double yhor_ME = y_ME; //Horizontal coordinates
    double zhor_ME = x_ME * cos(35*pi/180) + z_ME * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
    double az_ME  = atan2( yhor_ME, xhor_ME )*180/pi + 180; //local azimuth
    double alt_ME = atan2( zhor_ME, sqrt(xhor_ME*xhor_ME+yhor_ME*yhor_ME) )*180/pi; //local altitude
	double alt_topoc_ME = alt_ME - par_ME * cos(alt_ME*pi/180); //topocentric altitude
	double gclat_ME = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
    double rho_ME   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	double g_ME = atan( tan(gclat_ME*pi/180) / cos(HA_ME*pi/180) ); //an auxiliary angle
	double topRA_ME   = RA_ME  - par_ME * rho_ME * cos(gclat_ME*pi/180) * sin(HA_ME*pi/180) / cos(Dec_ME*pi/180);
    double topDecl_ME = Dec_ME - par_ME * rho_ME * sin(gclat_ME*pi/180) * sin(g_ME - Dec_ME*pi/180) / sin(g_ME);
	
    ////////////////////////////////Venus's Position////////////////////////////////
    double E_V = M_V + e_V*(180/pi) * sin(M_V*pi/180) * ( 1.0 + e_V * cos(M_V*pi/180) );  //Venus's eccentric anomaly E_S
    if(e_V>0.05){
        double E0_V=E_V;
        double E1_V=E_V + 2;
        while(fabs(E0_V-E1_V)>0.001){
            E0_V=E1_V;
            E1_V = E_V - ( E_V - e_V*(180/pi) * sin(E_V*pi/180) - M_V ) / ( 1 - e_V * cos(E_V*pi/180) );
        }
		E_V = E1_V;
    }
    double xv_V = a_V * (cos(E_V*pi/180) - e_V); //NOT IMPORTANT
    double yv_V = a_V * sqrt(1.0 - e_V*e_V) * sin(E_V*pi/180); //NOT IMPORTANT
    double v_V = atan2( yv_V, xv_V ); //Venus's true anomaly
    double r_V = sqrt( xv_V*xv_V + yv_V*yv_V ); //Venus's distance
    double xh_V = r_V * ( cos(N_V*pi/180) * cos(v_V+w_V*pi/180) - sin(N_V*pi/180) * sin(v_V+w_V*pi/180) * cos(i_V*pi/180) ); //Venus's heliocentric (Sun-centered) position in the ecliptic coordinate system
    double yh_V = r_V * ( sin(N_V*pi/180) * cos(v_V+w_V*pi/180) + cos(N_V*pi/180) * sin(v_V+w_V*pi/180) * cos(i_V*pi/180) ); //Venus's heliocentric (Sun-centered) position in the ecliptic coordinate system
    double zh_V = r_V * ( sin(v_V+w_V*pi/180) * sin(i_V*pi/180) ); //Venus's geocentric heliocentric (Sun-centered) position in the ecliptic coordinate system
	double lonecl_V = atan2( yh_V, xh_V )*180/pi; //Venus's ecliptic longitude
    double latecl_V = atan2( zh_V, sqrt(xh_V*xh_V+yh_V*yh_V) )*180/pi; //Venus's ecliptic latitude
	double xg_V = r_V * cos(lonecl_V*pi/180) * cos(latecl_V*pi/180) +xs; //geocentric position
	double yg_V = r_V * sin(lonecl_V*pi/180) * cos(latecl_V*pi/180) +ys; //geocentric position
	double zg_V = r_V * sin(latecl_V*pi/180); //geocentric position
	double xe_V = xg_V; //equatorial coordinates
    double ye_V = yg_V * cos(ecl*pi/180) - zg_V * sin(ecl*pi/180); //equatorial coordinates
    double ze_V = yg_V * sin(ecl*pi/180) + zg_V * cos(ecl*pi/180); //equatorial coordinates
    double RA_V  = atan2( ye_V, xe_V )*180/pi; //Right Ascension
    double Dec_V = atan2( ze_V, sqrt(xe_V*xe_V+ye_V*ye_V) )*180/pi; //Declination
	double par_V = (8.794/3600) / r_V;
	double HA_V = LST - RA_V; //hour angle of the Moon
	double x_V = cos(HA_V*pi/180) * cos(Dec_V*pi/180);
    double y_V = sin(HA_V*pi/180) * cos(Dec_V*pi/180);
    double z_V = sin(Dec_V*pi/180);
	double xhor_V = x_V * sin(35*pi/180) - z_V * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
    double yhor_V = y_V; //Horizontal coordinates
    double zhor_V = x_V * cos(35*pi/180) + z_V * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
    double az_V  = atan2( yhor_V, xhor_V )*180/pi + 180; //local azimuth
    double alt_V = atan2( zhor_V, sqrt(xhor_V*xhor_V+yhor_V*yhor_V) )*180/pi; //local altitude
	double alt_topoc_V = alt_V - par_V * cos(alt_V*pi/180); //topocentric altitude
	double gclat_V = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
    double rho_V   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	double g_V = atan( tan(gclat_V*pi/180) / cos(HA_V*pi/180) ); //an auxiliary angle
	double topRA_V   = RA_V  - par_V * rho_V * cos(gclat_V*pi/180) * sin(HA_V*pi/180) / cos(Dec_V*pi/180);
    double topDecl_V = Dec_V - par_V * rho_V * sin(gclat_V*pi/180) * sin(g_V - Dec_V*pi/180) / sin(g_V);
	
    ////////////////////////////////Mars's Position////////////////////////////////
    double E_MA = M_MA + e_MA*(180/pi) * sin(M_MA*pi/180) * ( 1.0 + e_MA * cos(M_MA*pi/180) );  //Mars's eccentric anomaly E_MA
    if(e_MA>0.05){
        double E0_MA=E_MA;
        double E1_MA=E_MA;
        while(fabs(E0_MA-E1_MA)>0.001){
            E0_MA=E1_MA;
            E1_MA = E_MA - ( E_MA - e_MA*(180/pi) * sin(E_MA*pi/180) - M_MA ) / ( 1 - e_MA * cos(E_MA*pi/180) );
        }
    }
    double xv_MA = a_MA * (cos(E_MA*pi/180) - e_MA); //NOT IMPORTANT
    double yv_MA = a_MA * sqrt(1.0 - e_MA*e_MA) * sin(E_MA*pi/180); //NOT IMPORTANT
    double v_MA = atan2( yv_MA, xv_MA ); //Mars's true anomaly
    double r_MA = sqrt( xv_MA*xv_MA + yv_MA*yv_MA ); //Mars's distance
    double xh_MA = r_MA * ( cos(N_MA*pi/180) * cos(v_MA+w_MA*pi/180) - sin(N_MA*pi/180) * sin(v_MA+w_MA*pi/180) * cos(i_MA*pi/180) ); //Mars's heliocentric (Sun-centered) position in the ecliptic coordinate system
    double yh_MA = r_MA * ( sin(N_MA*pi/180) * cos(v_MA+w_MA*pi/180) + cos(N_MA*pi/180) * sin(v_MA+w_MA*pi/180) * cos(i_MA*pi/180) ); //Mars's heliocentric (Sun-centered) position in the ecliptic coordinate system
    double zh_MA = r_MA * ( sin(v_MA+w_MA*pi/180) * sin(i_MA*pi/180) ); //Mars's heliocentric (Sun-centered) position in the ecliptic coordinate system
	double lonecl_MA = atan2( yh_MA, xh_MA )*180/pi; //Mars's ecliptic longitude
    double latecl_MA = atan2( zh_MA, sqrt(xh_MA*xh_MA+yh_MA*yh_MA) )*180/pi; //Mars's ecliptic latitude
	double xg_MA = r_MA * cos(lonecl_MA*pi/180) * cos(latecl_MA*pi/180) +xs; //geocentric position
	double yg_MA = r_MA * sin(lonecl_MA*pi/180) * cos(latecl_MA*pi/180) +ys; //geocentric position
	double zg_MA = r_MA * sin(latecl_MA*pi/180); //geocentric position
	double xe_MA = xg_MA; //equatorial coordinates
    double ye_MA = yg_MA * cos(ecl*pi/180) - zg_MA * sin(ecl*pi/180); //equatorial coordinates
    double ze_MA = yg_MA * sin(ecl*pi/180) + zg_MA * cos(ecl*pi/180); //equatorial coordinates
    double RA_MA  = atan2( ye_MA, xe_MA )*180/pi; //Right Ascension
    double Dec_MA = atan2( ze_MA, sqrt(xe_MA*xe_MA+ye_MA*ye_MA) )*180/pi; //Declination
	double par_MA = (8.794/3600) / r_MA;
	double HA_MA = LST - RA_MA; //hour angle of the Moon
	double x_MA = cos(HA_MA*pi/180) * cos(Dec_MA*pi/180);
    double y_MA = sin(HA_MA*pi/180) * cos(Dec_MA*pi/180);
    double z_MA = sin(Dec_MA*pi/180);
	double xhor_MA = x_MA * sin(35*pi/180) - z_MA * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
    double yhor_MA = y_MA; //Horizontal coordinates
    double zhor_MA = x_MA * cos(35*pi/180) + z_MA * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
    double az_MA  = atan2( yhor_MA, xhor_MA )*180/pi + 180; //local azimuth
    double alt_MA = atan2( zhor_MA, sqrt(xhor_MA*xhor_MA+yhor_MA*yhor_MA) )*180/pi; //local altitude
	double alt_topoc_MA = alt_MA - par_MA * cos(alt_MA*pi/180); //topocentric altitude
	double gclat_MA = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
    double rho_MA   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	double g_MA = atan( tan(gclat_MA*pi/180) / cos(HA_MA*pi/180) ); //an auxiliary angle
	double topRA_MA   = RA_MA  - par_MA * rho_MA * cos(gclat_MA*pi/180) * sin(HA_MA*pi/180) / cos(Dec_MA*pi/180);
    double topDecl_MA = Dec_MA - par_MA * rho_MA * sin(gclat_MA*pi/180) * sin(g_MA - Dec_MA*pi/180) / sin(g_MA);
	
    ////////////////////////////////Jupiter's Position////////////////////////////////
    double E_J = M_J + e_J*(180/pi) * sin(M_J*pi/180) * ( 1.0 + e_J * cos(M_J*pi/180) );  //Jupiter's eccentric anomaly E_J
    if(e_J>0.05){
        double E0_J=E_J;
        double E1_J=E_J;
        while(fabs(E0_J-E1_J)>0.001){
            E0_J=E1_J;
            E1_J = E_J - ( E_J - e_J*(180/pi) * sin(E_J*pi/180) - M_J ) / ( 1 - e_J * cos(E_J*pi/180) );
        }
    }
    double xv_J = a_J * (cos(E_J*pi/180) - e_J); //NOT IMPORTANT
    double yv_J = a_J * sqrt(1.0 - e_J*e_J) * sin(E_J*pi/180); //NOT IMPORTANT
    double v_J = atan2( yv_J, xv_J ); //Jupiter's true anomaly
    double r_J = sqrt( xv_J*xv_J + yv_J*yv_J ); //Mars's distance
    double xh_J = r_J * ( cos(N_J*pi/180) * cos(v_J+w_J*pi/180) - sin(N_J*pi/180) * sin(v_J+w_J*pi/180) * cos(i_J*pi/180) ); //Jupiter's heliocentric (Sun-centered) position in the ecliptic coordinate system
    double yh_J = r_J * ( sin(N_J*pi/180) * cos(v_J+w_J*pi/180) + cos(N_J*pi/180) * sin(v_J+w_J*pi/180) * cos(i_J*pi/180) ); //Jupiter's heliocentric (Sun-centered) position in the ecliptic coordinate system
    double zh_J = r_J * ( sin(v_J+w_J*pi/180) * sin(i_J*pi/180) ); //Jupiter's heliocentric (Sun-centered) position in the ecliptic coordinate system
	double lonecl_J = atan2( yh_J, xh_J )*180/pi; //Jupiter's ecliptic longitude
    double latecl_J = atan2( zh_J, sqrt(xh_J*xh_J+yh_J*yh_J) )*180/pi; //Jupiter's ecliptic latitude
	////////////////////////////////Saturn's Position////////////////////////////////
    double E_SA = M_SA + e_SA*(180/pi) * sin(M_SA*pi/180) * ( 1.0 + e_SA * cos(M_SA*pi/180) );  //Saturn's eccentric anomaly E_MA
    if(e_SA>0.05){
        double E0_SA=E_SA;
        double E1_SA=E_SA;
        while(fabs(E0_SA-E1_SA)>0.001){
            E0_SA=E1_SA;
            E1_SA = E_SA - ( E_SA - e_SA*(180/pi) * sin(E_SA*pi/180) - M_SA ) / ( 1 - e_SA * cos(E_SA*pi/180) );
        }
    }
    double xv_SA = a_SA * (cos(E_SA*pi/180) - e_SA); //NOT IMPORTANT
    double yv_SA = a_SA * sqrt(1.0 - e_MA*e_MA) * sin(E_MA*pi/180); //NOT IMPORTANT
    double v_SA = atan2( yv_MA, xv_MA ); //Saturn's true anomaly
    double r_SA = sqrt( xv_SA*xv_SA + yv_SA*yv_SA ); //Saturn's distance
    double xh_SA = r_SA * ( cos(N_SA*pi/180) * cos(v_SA+w_SA*pi/180) - sin(N_SA*pi/180) * sin(v_SA+w_SA*pi/180) * cos(i_SA*pi/180) ); //Saturn's heliocentric (Sun-centered) position in the ecliptic coordinate system
    double yh_SA = r_SA * ( sin(N_SA*pi/180) * cos(v_SA+w_SA*pi/180) + cos(N_SA*pi/180) * sin(v_SA+w_SA*pi/180) * cos(i_SA*pi/180) ); //Saturn's heliocentric (Sun-centered) position in the ecliptic coordinate system
    double zh_SA = r_SA * ( sin(v_SA+w_SA*pi/180) * sin(i_SA*pi/180) ); //Saturn's heliocentric (Sun-centered) position in the ecliptic coordinate system
	double lonecl_SA = atan2( yh_SA, xh_SA )*180/pi; //Saturn's ecliptic longitude
    double latecl_SA = atan2( zh_SA, sqrt(xh_SA*xh_SA+yh_SA*yh_SA) )*180/pi; //Saturn's ecliptic latitude
	////////////////////////////////Uranus's Position////////////////////////////////
    double E_U = M_U + e_U*(180/pi) * sin(M_U*pi/180) * ( 1.0 + e_U * cos(M_U*pi/180) );  //Uranus's eccentric anomaly E_MA
    if(e_U>0.05){
        double E0_U=E_U;
        double E1_U=E_U;
        while(fabs(E0_U-E1_U)>0.001){
            E0_U=E1_U;
            E1_U = E_U - ( E_U - e_U*(180/pi) * sin(E_U*pi/180) - M_U ) / ( 1 - e_U * cos(E_U*pi/180) );
        }
    }
    double xv_U = a_U * (cos(E_U*pi/180) - e_U); //NOT IMPORTANT
    double yv_U = a_U * sqrt(1.0 - e_U*e_U) * sin(E_MA*pi/180); //NOT IMPORTANT
    double v_U = atan2( yv_U, xv_U ); //Uranus's true anomaly
    double r_U = sqrt( xv_U*xv_U + yv_U*yv_U ); //Uranus's distance
    double xh_U = r_U * ( cos(N_U*pi/180) * cos(v_U+w_U*pi/180) - sin(N_U*pi/180) * sin(v_U+w_U*pi/180) * cos(i_U*pi/180) ); //Uranus's heliocentric (Sun-centered) position in the ecliptic coordinate system
    double yh_U = r_U * ( sin(N_U*pi/180) * cos(v_U+w_U*pi/180) + cos(N_U*pi/180) * sin(v_U+w_U*pi/180) * cos(i_U*pi/180) ); //Uranus's heliocentric (Sun-centered) position in the ecliptic coordinate system
    double zh_U = r_U * ( sin(v_U+w_U*pi/180) * sin(i_U*pi/180) ); //Uranus's heliocentric (Sun-centered) position in the ecliptic coordinate system
	double lonecl_U = atan2( yh_U, xh_U )*180/pi; //Uranus's ecliptic longitude
    double latecl_U = atan2( zh_U, sqrt(xh_U*xh_U+yh_U*yh_U) )*180/pi; //Uranus's ecliptic latitude
	////////////////////////////////Neptune's Position////////////////////////////////
    double E_N = M_N + e_N*(180/pi) * sin(M_N*pi/180) * ( 1.0 + e_N * cos(M_N*pi/180) );  //Neptune's eccentric anomaly E_MA
    if(e_N>0.05){
        double E0_N=E_N;
        double E1_N=E_N;
        while(fabs(E0_N-E1_N)>0.001){
            E0_N=E1_N;
            E1_N = E_N - ( E_N - e_N*(180/pi) * sin(E_N*pi/180) - M_N ) / ( 1 - e_N * cos(E_N*pi/180) );
        }
    }
    double xv_N = a_N * (cos(E_N*pi/180) - e_N); //NOT IMPORTANT
    double yv_N = a_N * sqrt(1.0 - e_N*e_N) * sin(E_N*pi/180); //NOT IMPORTANT
    double v_N = atan2( yv_N, xv_N ); //Neptune's true anomaly
    double r_N = sqrt( xv_N*xv_N + yv_N*yv_N ); //Neptune's distance
    double xh_N = r_N * ( cos(N_N*pi/180) * cos(v_N+w_N*pi/180) - sin(N_N*pi/180) * sin(v_N+w_N*pi/180) * cos(i_N*pi/180) ); //Neptune's heliocentric (Sun-centered) position in the ecliptic coordinate system
    double yh_N = r_N * ( sin(N_N*pi/180) * cos(v_N+w_N*pi/180) + cos(N_N*pi/180) * sin(v_N+w_N*pi/180) * cos(i_N*pi/180) ); //Neptune's heliocentric (Sun-centered) position in the ecliptic coordinate system
    double zh_N = r_N * ( sin(v_N+w_N*pi/180) * sin(i_N*pi/180) ); //Neptune's heliocentric (Sun-centered) position in the ecliptic coordinate system
	double lonecl_N = atan2( yh_N, xh_N )*180/pi; //Neptune's ecliptic longitude
    double latecl_N = atan2( zh_N, sqrt(xh_N*xh_N+yh_N*yh_N) )*180/pi; //Neptune's ecliptic latitude		double xg_N = xh_N + xs; //geocentric position
	double xg_N = r_N * cos(lonecl_N*pi/180) * cos(latecl_N*pi/180) +xs; //geocentric position
	double yg_N = r_N * sin(lonecl_N*pi/180) * cos(latecl_N*pi/180) +ys; //geocentric position
	double zg_N = r_N * sin(latecl_N*pi/180); //geocentric position
	double xe_N = xg_N; //equatorial coordinates
    double ye_N = yg_N * cos(ecl*pi/180) - zg_N * sin(ecl*pi/180); //equatorial coordinates
    double ze_N = yg_N * sin(ecl*pi/180) + zg_N * cos(ecl*pi/180); //equatorial coordinates
    double RA_N  = atan2( ye_N, xe_N )*180/pi; //Right Ascension
    double Dec_N = atan2( ze_N, sqrt(xe_N*xe_N+ye_N*ye_N) )*180/pi; //Declination
	double par_N = (8.794/3600) / r_N;
	double HA_N = LST - RA_N; //hour angle of the Moon
	double x_N = cos(HA_N*pi/180) * cos(Dec_N*pi/180);
    double y_N = sin(HA_N*pi/180) * cos(Dec_N*pi/180);
    double z_N = sin(Dec_N*pi/180);
	double xhor_N = x_N * sin(35*pi/180) - z_N * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
    double yhor_N = y_N; //Horizontal coordinates
    double zhor_N = x_N * cos(35*pi/180) + z_N * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
    double az_N  = atan2( yhor_N, xhor_N )*180/pi + 180; //local azimuth
    double alt_N = atan2( zhor_N, sqrt(xhor_N*xhor_N+yhor_N*yhor_N) )*180/pi; //local altitude
	double (8.794/3600) / r_N;
	double alt_topoc_N = alt_N - par_N * cos(alt_N*pi/180); //topocentric altitude
	double gclat_N = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
    double rho_N   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	double g_N = atan( tan(gclat_N*pi/180) / cos(HA_N*pi/180) ); //an auxiliary angle
	double topRA_N   = RA_N  - par_N * rho_N * cos(gclat_N*pi/180) * sin(HA_N*pi/180) / cos(Dec_N*pi/180);
    double topDecl_N = Dec_N - par_N * rho_N * sin(gclat_N*pi/180) * sin(g_N - Dec_N*pi/180) / sin(g_N);
	
	////////////////////////////////Pluto's Position////////////////////////////////
    double S  =   50.03  +  0.033459652 * d;
    double P  =  238.95  +  0.003968789 * d;
	double lonecl_P = 238.9508  +  0.00400703 * d
            - 19.799 * sin(P*pi/180)     + 19.848 * cos(P*pi/180)
             + 0.897 * sin(2*P*pi/180)    - 4.956 * cos(2*P*pi/180)
             + 0.610 * sin(3*P*pi/180)    + 1.211 * cos(3*P*pi/180)
             - 0.341 * sin(4*P*pi/180)    - 0.190 * cos(4*P*pi/180)
             + 0.128 * sin(5*P*pi/180)    - 0.034 * cos(5*P*pi/180)
             - 0.038 * sin(6*P*pi/180)    + 0.031 * cos(6*P*pi/180)
             + 0.020 * sin(S*pi/180-P*pi/180)    - 0.010 * cos(S*pi/180-P*pi/180);
	 double latecl_P = 453 * sin(P*pi/180)     - 14.975 * cos(P*pi/180)
             + 3.527 * sin(2*P*pi/180)    + 1.673 * cos(2*P*pi/180)
             - 1.051 * sin(3*P*pi/180)    + 0.328 * cos(3*P*pi/180)
             + 0.179 * sin(4*P*pi/180)    - 0.292 * cos(4*P*pi/180)
             + 0.019 * sin(5*P*pi/180)    + 0.100 * cos(5*P*pi/180)
             - 0.031 * sin(6*P*pi/180)    - 0.026 * cos(6*P*pi/180)
                                   + 0.011 * cos(S*pi/180-P*pi/180);
	 double r_P     =  40.72
           + 6.68 * sin(P*pi/180)       + 6.90 * cos(P*pi/180)
           - 1.18 * sin(2*P*pi/180)     - 0.03 * cos(2*P*pi/180)
           + 0.15 * sin(3*P*pi/180)     - 0.14 * cos(3*P*pi/180);
	 
	while(lonecl_P>90)
		lonecl_P = lonecl_P - 90;
	while(lonecl_P<-90)
		lonecl_P = lonecl_P + 90;

	while(latecl_P>90)
		latecl_P = latecl_P - 90;
	while(latecl_P<-90)
		latecl_P = latecl_P + 90;
    double xg_P = r_P * cos(lonecl_P*pi/180) * cos(latecl_P*pi/180) +xs; //geocentric position
	double yg_P = r_P * sin(lonecl_P*pi/180) * cos(latecl_P*pi/180) +ys; //geocentric position
	double zg_P = r_P * sin(latecl_P*pi/180); //geocentric position
	double xe_P = xg_P; //equatorial coordinates
    double ye_P = yg_P * cos(ecl*pi/180) - zg_P * sin(ecl*pi/180); //equatorial coordinates
    double ze_P = yg_P * sin(ecl*pi/180) + zg_P * cos(ecl*pi/180); //equatorial coordinates
    double RA_P  = atan2( ye_P, xe_P )*180/pi; //Right Ascension
    double Dec_P = atan2( ze_P, sqrt(xe_P*xe_P+ye_P*ye_P) )*180/pi; //Declination
	double par_P = (8.794/3600) / r_P;
	double HA_P = LST - RA_P; //hour angle of the Moon
	double x_P = cos(HA_P*pi/180) * cos(Dec_P*pi/180);
    double y_P = sin(HA_P*pi/180) * cos(Dec_P*pi/180);
    double z_P = sin(Dec_P*pi/180);
	double xhor_P = x_P * sin(35*pi/180) - z_P * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
    double yhor_P = y_P; //Horizontal coordinates
    double zhor_P = x_P * cos(35*pi/180) + z_P * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
    double az_P  = atan2( yhor_P, xhor_P )*180/pi + 180; //local azimuth
    double alt_P = atan2( zhor_P, sqrt(xhor_P*xhor_P+yhor_P*yhor_P) )*180/pi; //local altitude
	double alt_topoc_P = alt_P - par_P * cos(alt_P*pi/180); //topocentric altitude
	double gclat_P = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
    double rho_P   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	double g_P = atan( tan(gclat_P*pi/180) / cos(HA_P*pi/180) ); //an auxiliary angle
	double topRA_P   = RA_P  - par_P * rho_P * cos(gclat_P*pi/180) * sin(HA_P*pi/180) / cos(Dec_P*pi/180);
    double topDecl_P = Dec_P - par_P * rho_P * sin(gclat_P*pi/180) * sin(g_P - Dec_P*pi/180) / sin(g_P);
	
    ////////////////////////////////////////////////////////////////

	////////////////////////////////Perturbations of the Moon////////////////////////////////
	lonecl_M = lonecl_M -1.274 * sin(M_M*pi/180 - 2*D*pi/180)          
    +0.658 * sin(2*D*pi/180)               
    -0.186 * sin(M_S*pi/180)                
    -0.059 * sin(2*M_M*pi/180 - 2*D*pi/180)
    -0.057 * sin(M_M*pi/180 - 2*D*pi/180 + M_S*pi/180)
    +0.053 * sin(M_M*pi/180 + 2*D*pi/180)
    +0.046 * sin(2*D*pi/180 - M_S*pi/180)
    +0.041 * sin(M_M*pi/180 - M_S*pi/180)
    -0.035 * sin(D*pi/180)                 
    -0.031 * sin(M_M*pi/180 + M_S*pi/180)
    -0.015 * sin(2*F*pi/180 - 2*D*pi/180)
    +0.011 * sin(M_M*pi/180 - 4*D*pi/180);

	latecl_M = latecl_M -0.173 * sin(F*pi/180 - 2*D*pi/180)
    -0.055 * sin(M_M*pi/180 - F*pi/180 - 2*D*pi/180)
    -0.046 * sin(M_M*pi/180 + F*pi/180 - 2*D*pi/180)
    +0.033 * sin(F*pi/180 + 2*D*pi/180)
    +0.017 * sin(2*M_M*pi/180 + F*pi/180);

	r_M = r_M -0.58 * cos(M_M*pi/180 - 2*D*pi/180)
    -0.46 * cos(2*D*pi/180);

	double xg_M = r_M * cos(lonecl_M*pi/180) * cos(latecl_M*pi/180); //xh_M + xs; //geocentric position
	double yg_M = r_M * sin(lonecl_M*pi/180) * cos(latecl_M*pi/180); //yh_M + ys; //geocentric position
	double zg_M = r_M * sin(latecl_M*pi/180); //zh_M; //geocentric position
	double xe_M = xg_M; //equatorial coordinates
    double ye_M = yg_M * cos(ecl*pi/180) - zg_M * sin(ecl*pi/180); //equatorial coordinates
    double ze_M = yg_M * sin(ecl*pi/180) + zg_M * cos(ecl*pi/180); //equatorial coordinates
    double RA_M  = atan2( ye_M, xe_M )*180/pi; //Right Ascension
    double Dec_M = atan2( ze_M, sqrt(xe_M*xe_M+ye_M*ye_M) )*180/pi; //Declination
	double HA_M = LST - RA_M; //hour angle of the Moon
	double x_M = cos(HA_M*pi/180) * cos(Dec_M*pi/180);
    double y_M = sin(HA_M*pi/180) * cos(Dec_M*pi/180);
    double z_M = sin(Dec_M*pi/180);
	double xhor_M = x_M * sin(35*pi/180) - z_M * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
    double yhor_M = y_M; //Horizontal coordinates
    double zhor_M = x_M * cos(35*pi/180) + z_M * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
    double az_M  = atan2( yhor_M, xhor_M )*180/pi + 180; //local azimuth
    double alt_M = atan2( zhor_M, sqrt(xhor_M*xhor_M+yhor_M*yhor_M) )*180/pi; //local altitude
	double mpar = asin( 1/r_M )*180/pi; //Moon's Paralax
	double alt_topoc_M = alt_M - mpar * cos(alt_M*pi/180); //topocentric altitude
	double gclat_M = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
    double rho_M   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	//double gclat = 35;
	//double rho = 1.0;
	double g_M = atan( tan(gclat_M*pi/180) / cos(HA_M*pi/180) ); //an auxiliary angle
	double topRA_M   = RA_M  - mpar * rho_M * cos(gclat_M*pi/180) * sin(HA_M*pi/180) / cos(Dec_M*pi/180);
    double topDecl_M = Dec_M - mpar * rho_M * sin(gclat_M*pi/180) * sin(g_M - Dec_M*pi/180) / sin(g_M);
	////////////////////////////////////////////////////////////////

	////////////////////////////////Perturbations of Jupiter, Saturn and Uranus////////////////////////////////
	lonecl_J = lonecl_J -0.332 * sin(2*M_J*pi/180 - 5*M_SA*pi/180 - 67.6*pi/180 )
    -0.056 * sin(2*M_J*pi/180 - 2*M_SA*pi/180 + 21*pi/180 )
    +0.042 * sin(3*M_J*pi/180 - 5*M_SA*pi/180 + 21*pi/180 )
    -0.036 * sin(M_J*pi/180 - 2*M_SA*pi/180)
    +0.022 * cos(M_J*pi/180 - M_SA*pi/180)
    +0.023 * sin(2*M_J*pi/180 - 3*M_SA + 52*pi/180 )
    -0.016 * sin(M_J*pi/180 - 5*M_SA*pi/180 - 69*pi/180 );

	lonecl_SA = lonecl_SA +0.812 * sin(2*M_J*pi/180 - 5*M_SA*pi/180 - 67.6*pi/180 )
    -0.229 * cos(2*M_J*pi/180 - 4*M_SA*pi/180 - 2*pi/180 )
    +0.119 * sin(M_J*pi/180 - 2*M_SA*pi/180 - 3*pi/180 )
    +0.046 * sin(2*M_J*pi/180 - 6*M_SA*pi/180 - 69*pi/180 )
    +0.014 * sin(M_J*pi/180 - 3*M_SA*pi/180 + 32*pi/180 );

	latecl_SA = latecl_SA -0.020 * cos(2*M_J*pi/180 - 4*M_SA*pi/180 - 2*pi/180 )
    +0.018 * sin(2*M_J*pi/180 - 6*M_SA*pi/180 - 49*pi/180 );

	lonecl_U = lonecl_U + 0.040 * sin(M_SA*pi/180 - 2*M_U*pi/180 + 6*pi/180)
    +0.035 * sin(M_SA*pi/180 - 3*M_U*pi/180 + 33*pi/180)
    -0.015 * sin(M_J*pi/180 - M_U*pi/180 + 20*pi/180);

	double xg_J = r_J * cos(lonecl_J*pi/180) * cos(latecl_J*pi/180) +xs; //geocentric position
	double yg_J = r_J * sin(lonecl_J*pi/180) * cos(latecl_J*pi/180) +ys; //geocentric position
	double zg_J = r_J * sin(latecl_J*pi/180); //geocentric position
	double xe_J = xg_J; //equatorial coordinates
    double ye_J = yg_J * cos(ecl*pi/180) - zg_J * sin(ecl*pi/180); //equatorial coordinates
    double ze_J = yg_J * sin(ecl*pi/180) + zg_J * cos(ecl*pi/180); //equatorial coordinates
    double RA_J  = atan2( ye_J, xe_J )*180/pi; //Right Ascension
    double Dec_J = atan2( ze_J, sqrt(xe_J*xe_J+ye_J*ye_J) )*180/pi; //Declination
	double par_J = (8.794/3600) / r_J;
	double HA_J = LST - RA_J; //hour angle of the Moon
	double x_J = cos(HA_J*pi/180) * cos(Dec_J*pi/180);
    double y_J = sin(HA_J*pi/180) * cos(Dec_J*pi/180);
    double z_J = sin(Dec_J*pi/180);
	double xhor_J = x_J * sin(35*pi/180) - z_J * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
    double yhor_J = y_J; //Horizontal coordinates
    double zhor_J = x_J * cos(35*pi/180) + z_J * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
    double az_J  = atan2( yhor_J, xhor_J )*180/pi + 180; //local azimuth
    double alt_J = atan2( zhor_J, sqrt(xhor_J*xhor_J+yhor_J*yhor_J) )*180/pi; //local altitude
	double alt_topoc_J = alt_J - par_J * cos(alt_J*pi/180); //topocentric altitude
	double gclat_J = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
    double rho_J   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	double g_J = atan( tan(gclat_J*pi/180) / cos(HA_J*pi/180) ); //an auxiliary angle
	double topRA_J   = RA_J  - par_J * rho_J * cos(gclat_J*pi/180) * sin(HA_J*pi/180) / cos(Dec_J*pi/180);
    double topDecl_J = Dec_J - par_J * rho_J * sin(gclat_J*pi/180) * sin(g_J - Dec_J*pi/180) / sin(g_J);
	
	
	double xg_SA = r_SA * cos(lonecl_SA*pi/180) * cos(latecl_SA*pi/180) +xs; //geocentric position
	double yg_SA = r_SA * sin(lonecl_SA*pi/180) * cos(latecl_SA*pi/180) +ys; //geocentric position
	double zg_SA = r_SA * sin(latecl_SA*pi/180); //geocentric position
	double xe_SA = xg_SA; //equatorial coordinates
    double ye_SA = yg_SA * cos(ecl*pi/180) - zg_SA * sin(ecl*pi/180); //equatorial coordinates
    double ze_SA = yg_SA * sin(ecl*pi/180) + zg_SA * cos(ecl*pi/180); //equatorial coordinates
    double RA_SA  = atan2( ye_SA, xe_SA )*180/pi; //Right Ascension
    double Dec_SA = atan2( ze_SA, sqrt(xe_SA*xe_SA+ye_SA*ye_SA) )*180/pi; //Declination
    double par_SA = (8.794/3600) / r_SA;
	double HA_SA = LST - RA_SA; //hour angle of the Moon
	double x_SA = cos(HA_SA*pi/180) * cos(Dec_SA*pi/180);
    double y_SA = sin(HA_SA*pi/180) * cos(Dec_SA*pi/180);
    double z_SA = sin(Dec_SA*pi/180);
	double xhor_SA = x_SA * sin(35*pi/180) - z_SA * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
    double yhor_SA = y_SA; //Horizontal coordinates
    double zhor_SA = x_SA * cos(35*pi/180) + z_SA * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
    double az_SA  = atan2( yhor_SA, xhor_SA )*180/pi + 180; //local azimuth
    double alt_SA = atan2( zhor_SA, sqrt(xhor_SA*xhor_SA+yhor_SA*yhor_SA) )*180/pi; //local altitude
	double alt_topoc_SA = alt_SA - par_SA * cos(alt_SA*pi/180); //topocentric altitude
	double gclat_SA = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
    double rho_SA   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	double g_SA = atan( tan(gclat_SA*pi/180) / cos(HA_SA*pi/180) ); //an auxiliary angle
	double topRA_SA   = RA_SA  - par_SA * rho_SA * cos(gclat_SA*pi/180) * sin(HA_SA*pi/180) / cos(Dec_SA*pi/180);
    double topDecl_SA = Dec_SA - par_SA * rho_SA * sin(gclat_SA*pi/180) * sin(g_SA - Dec_SA*pi/180) / sin(g_SA);
	

	double xg_U = r_U * cos(lonecl_U*pi/180) * cos(latecl_U*pi/180) +xs; //geocentric position
	double yg_U = r_U * sin(lonecl_U*pi/180) * cos(latecl_U*pi/180) +ys; //geocentric position
	double zg_U = r_U * sin(latecl_U*pi/180); //geocentric position
	double xe_U = xg_U; //equatorial coordinates
    double ye_U = yg_U * cos(ecl*pi/180) - zg_U * sin(ecl*pi/180); //equatorial coordinates
    double ze_U = yg_U * sin(ecl*pi/180) + zg_U * cos(ecl*pi/180); //equatorial coordinates
    double RA_U  = atan2( ye_U, xe_U )*180/pi; //Right Ascension
    double Dec_U = atan2( ze_U, sqrt(xe_U*xe_U+ye_U*ye_U) )*180/pi; //Declination
    double par_U = (8.794/3600) / r_U;
	double HA_U = LST - RA_U; //hour angle of the Moon
	double x_U = cos(HA_U*pi/180) * cos(Dec_U*pi/180);
    double y_U = sin(HA_U*pi/180) * cos(Dec_U*pi/180);
    double z_U = sin(Dec_U*pi/180);
	double xhor_U = x_U * sin(35*pi/180) - z_U * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
    double yhor_U = y_U; //Horizontal coordinates
    double zhor_U = x_U * cos(35*pi/180) + z_U * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
    double az_U  = atan2( yhor_U, xhor_U )*180/pi + 180; //local azimuth
    double alt_U = atan2( zhor_U, sqrt(xhor_U*xhor_U+yhor_U*yhor_U) )*180/pi; //local altitude
	double alt_topoc_U = alt_U - par_U * cos(alt_U*pi/180); //topocentric altitude
	double gclat_U = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
    double rho_U   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	double g_U = atan( tan(gclat_U*pi/180) / cos(HA_U*pi/180) ); //an auxiliary angle
	double topRA_U   = RA_U  - par_U * rho_U * cos(gclat_U*pi/180) * sin(HA_U*pi/180) / cos(Dec_U*pi/180);
    double topDecl_U = Dec_U - par_U * rho_U * sin(gclat_U*pi/180) * sin(g_U - Dec_U*pi/180) / sin(g_U);
	
    ////////////////////////////////////////////////////////////////


    ////////////////////////////////Precession////////////////////////////////
	//lon_corr = 3.82394E-5 * ( 365.2422 * ( Epoch - 2000.0 ) - d ) //uncoM_Ment this if you need to do your calculations for some standard epoch, such as 1950.0 or 2000.0. 
    ////////////////////////////////////////////////////////////////
	while (RA_S>360)
		RA_S = RA_S -360;
	while (RA_S<0)
		RA_S = RA_S +360;
	while (Dec_S<-90)
		Dec_S = Dec_S + 90;
	while (Dec_S>90)
		Dec_S = Dec_S - 90;
		while (RA_M>360)
		RA_M = RA_M -360;
	while (RA_M<0)
		RA_M = RA_M +360;
	while (Dec_M<-90)
		Dec_M = Dec_M + 90;
	while (Dec_M>90)
		Dec_M = Dec_M - 90;
	while (RA_ME>360)
		RA_ME = RA_ME -360;
	while (RA_ME<0)
		RA_ME = RA_ME +360;
	while (Dec_ME<-90)
		Dec_ME = Dec_ME + 90;
	while (Dec_ME>90)
		Dec_ME = Dec_ME- 90;
		while (RA_V>360)
		RA_V = RA_V -360;
	while (RA_V<0)
		RA_V = RA_V +360;
	while (Dec_V<-90)
		Dec_V = Dec_V + 90;
	while (Dec_V>90)
		Dec_V = Dec_V - 90;
		while (RA_MA>360)
		RA_MA = RA_MA -360;
	while (RA_MA<0)
		RA_MA = RA_MA +360;
	while (Dec_MA<-90)
		Dec_MA = Dec_MA + 90;
	while (Dec_MA>90)
		Dec_MA = Dec_MA - 90;
		while (RA_J>360)
		RA_J = RA_J -360;
	while (RA_J<0)
		RA_J = RA_J +360;
	while (Dec_J<-90)
		Dec_J = Dec_J + 90;
	while (Dec_J>90)
		Dec_J = Dec_J - 90;
		while (RA_SA>360)
		RA_SA = RA_SA -360;
	while (RA_SA<0)
		RA_SA = RA_SA +360;
	while (Dec_SA<-90)
		Dec_SA = Dec_SA + 90;
	while (Dec_SA>90)
		Dec_SA = Dec_SA - 90;
		while (RA_U>360)
		RA_U = RA_U -360;
	while (RA_U<0)
		RA_U = RA_U +360;
	while (Dec_U<-90)
		Dec_U = Dec_U + 90;
	while (Dec_U>90)
		Dec_U = Dec_U - 90;
		while (RA_N>360)
		RA_N = RA_N -360;
	while (RA_N<0)
		RA_N = RA_N +360;
	while (Dec_N<-90)
		Dec_N = Dec_N + 90;
	while (Dec_N>90)
		Dec_N = Dec_N - 90;
		while (RA_P>360)
		RA_P = RA_P -360;
	while (RA_P<0)
		RA_P = RA_P +360;
	while (Dec_P<-90)
		Dec_P = Dec_P + 90;
	while (Dec_P>90)
		Dec_P = Dec_P - 90;
			
		while (topRA_M>360)
		topRA_M = topRA_M -360;
	while (topRA_M<0)
		topRA_M = topRA_M +360;
	while (topDecl_M<-90)
		topDecl_M = topDecl_M + 90;
	while (topDecl_M>90)
		topDecl_M = topDecl_M - 90;
	while (topRA_ME>360)
		topRA_ME = topRA_ME -360;
	while (topRA_ME<0)
		topRA_ME = topRA_ME +360;
	while (topDecl_ME<-90)
		topDecl_ME = topDecl_ME + 90;
	while (topDecl_ME>90)
		topDecl_ME = topDecl_ME - 90;
		
		while (topRA_V>360)
		topRA_V = topRA_V -360;
	while (topRA_V<0)
		topRA_V = topRA_V +360;
	while (topDecl_V<-90)
		topDecl_V = topDecl_V + 90;
	while (topDecl_V>90)
		topDecl_V = topDecl_V - 90;
		while (topRA_MA>360)
		topRA_MA = topRA_MA -360;
	while (topRA_MA<0)
		topRA_MA = topRA_MA +360;
	while (topDecl_MA<-90)
		topDecl_MA = topDecl_MA + 90;
	while (topDecl_MA>90)
		topDecl_MA = topDecl_MA - 90;
		while (topRA_J>360)
		topRA_J = topRA_J -360;
	while (topRA_J<0)
		topRA_J = topRA_J +360;
	while (topDecl_J<-90)
		topDecl_J = topDecl_J + 90;
	while (topDecl_J>90)
		topDecl_J = topDecl_J - 90;
		while (topRA_SA>360)
		topRA_SA = topRA_SA -360;
	while (topRA_SA<0)
		topRA_SA = topRA_SA +360;
	while (topDecl_SA<-90)
		topDecl_SA = topDecl_SA + 90;
	while (topDecl_SA>90)
		topDecl_SA = topDecl_SA - 90;
		while (topRA_U>360)
		topRA_U = topRA_U -360;
	while (topRA_U<0)
		topRA_U = topRA_U +360;
	while (topDecl_U<-90)
		topDecl_U = topDecl_U + 90;
	while (topDecl_U>90)
		topDecl_U = topDecl_U - 90;
	while (topRA_N>360)
		topRA_N = topRA_N -360;
	while (topRA_N<0)
		topRA_N = topRA_N +360;
	while (topDecl_N<-90)
		topDecl_N = topDecl_N + 90;
	while (topDecl_N>90)
		topDecl_N = topDecl_N - 90;
		while (topRA_P>360)
		topRA_P = topRA_P -360;
	while (topRA_P<0)
		topRA_P = topRA_P +360;
	while (topDecl_P<-90)
		topDecl_P = topDecl_P + 90;
	while (topDecl_P>90)
		topDecl_P = topDecl_P - 90;
		
	cout<<"For the Sun:"<<endl<<"*************************************************************************"<<endl;
	cout<<"The eccentricity (e): "<<e_S<<endl;
    cout<<"The mean distance from the earth (a): "<<a_S<<endl;
	cout<<"The inclination to the ecliptic (i): "<<i_S<<endl;
	cout<<"The longitude of the ascending node (Omega): "<<N_S<<endl;
	cout<<"The argument of perihelion (w): "<<w_S<<endl;
	cout<<"The mean anomaly (M): "<<M_S<<endl;
	cout<<"The obliquity of the ecliptic (epsilon): "<<ecl<<endl;
	cout<<"The eccentric anomaly (E): "<<E_S<<endl;
	cout<<"The true longitude (lambda): "<<lonsun<<endl;
	cout<<"The true latitude (beta): "<<0.0<<endl;
	cout<<"The distance (r): "<<rs<<endl;
	cout<<"Geocentric coordinates:"<<endl<<"x = "<<xs<<endl<<"y = "<<ys<<endl<<"z = "<<0.0<<endl;
	cout<<"Equatorial coordinates:"<<endl<<"x = "<<xes<<endl<<"y = "<<yes<<endl<<"z = "<<zes<<endl<<"alpha = "<<RA_S<<endl<<"delta = "<<Dec_S<<endl;
	cout<<"*************************************************************************";

	cout<<endl<<"For the Moon:"<<endl<<"*************************************************************************"<<endl;
	cout<<"The eccentricity (e): "<<e_M<<endl;
    cout<<"The mean distance from the earth (a): "<<a_M<<" x R(earth)"<<endl;
	cout<<"The inclination to the ecliptic (i): "<<i_M<<endl;
	cout<<"The longitude of the ascending node (Omega): "<<N_M<<endl;
	cout<<"The argument of perihelion (w): "<<w_M<<endl;
	cout<<"The mean anomaly (M): "<<M_M<<endl;
	cout<<"The eccentric anomaly (E): "<<E_M<<endl;
	cout<<"The true longitude (lambda): "<<lonecl_M<<endl;
	cout<<"The true latitude (beta): "<<latecl_M<<endl;
	cout<<"The distance (r): "<<r_M<<" x R(earth)"<<endl;
	cout<<"Geocentric coordinates:"<<endl<<"x = "<<xg_M<<" x R(earth)"<<endl<<"y = "<<yg_M<<" x R(earth)"<<endl<<"z = "<<zg_M<<" x R(earth)"<<endl;
	cout<<"Equatorial coordinates:"<<endl<<"x = "<<xe_M<<" x R(earth)"<<endl<<"y = "<<ye_M<<" x R(earth)"<<endl<<"z = "<<ze_M<<" x R(earth)"<<endl<<"alpha = "<<RA_M<<endl<<"delta = "<<Dec_M<<endl;
	cout<<"Topocentric coordinates"<<endl<<"alpha = "<<topRA_M<<endl<<"delta = "<<topDecl_M<<endl;
	cout<<"*************************************************************************";

	cout<<endl<<"For Mercury:"<<endl<<"*************************************************************************"<<endl;
	cout<<"The eccentricity (e): "<<e_ME<<endl;
    cout<<"The mean distance from the earth (a): "<<a_ME<<endl;
	cout<<"The inclination to the ecliptic (i): "<<i_ME<<endl;
	cout<<"The longitude of the ascending node (Omega): "<<N_ME<<endl;
	cout<<"The argument of perihelion (w): "<<w_ME<<endl;
	cout<<"The mean anomaly (M): "<<M_ME<<endl;
	cout<<"The eccentric anomaly (E): "<<E_ME<<endl;
	cout<<"The true longitude (lambda): "<<lonecl_ME<<endl;
	cout<<"The true latitude (beta): "<<latecl_ME<<endl;
	cout<<"The distance (r): "<<r_ME<<endl;
	cout<<"Geocentric coordinates:"<<endl<<"x = "<<xg_ME<<endl<<"y = "<<yg_ME<<endl<<"z = "<<zg_ME<<endl;
	cout<<"Equatorial coordinates:"<<endl<<"x = "<<xe_ME<<endl<<"y = "<<ye_ME<<endl<<"z = "<<ze_ME<<endl<<"alpha = "<<RA_ME<<endl<<"delta = "<<Dec_ME<<endl;
	cout<<"Topocentric coordinates"<<endl<<"alpha = "<<topRA_ME<<endl<<"delta = "<<topDecl_ME<<endl;
	cout<<"*************************************************************************";

	cout<<endl<<"For Venus:"<<endl<<"*************************************************************************"<<endl;
	cout<<"The eccentricity (e): "<<e_V<<endl;
    cout<<"The mean distance from the earth (a): "<<a_V<<endl;
	cout<<"The inclination to the ecliptic (i): "<<i_V<<endl;
	cout<<"The longitude of the ascending node (Omega): "<<N_V<<endl;
	cout<<"The argument of perihelion (w): "<<w_V<<endl;
	cout<<"The mean anomaly (M): "<<M_V<<endl;
	cout<<"The eccentric anomaly (E): "<<E_V<<endl;
	cout<<"The true longitude (lambda): "<<lonecl_V<<endl;
	cout<<"The true latitude (beta): "<<latecl_V<<endl;
	cout<<"The distance (r): "<<r_V<<endl;
	cout<<"Geocentric coordinates:"<<endl<<"x = "<<xg_V<<endl<<"y = "<<yg_V<<endl<<"z = "<<zg_V<<endl;
	cout<<"Equatorial coordinates:"<<endl<<"x = "<<xe_V<<endl<<"y = "<<ye_V<<endl<<"z = "<<ze_V<<endl<<"alpha = "<<RA_V<<endl<<"delta = "<<Dec_V<<endl;
	cout<<"Topocentric coordinates"<<endl<<"alpha = "<<topRA_V<<endl<<"delta = "<<topDecl_V<<endl;
	
	cout<<"*************************************************************************";

	cout<<endl<<"For Mars:"<<endl<<"*************************************************************************"<<endl;
	cout<<"The eccentricity (e): "<<e_MA<<endl;
    cout<<"The mean distance from the earth (a): "<<a_MA<<endl;
	cout<<"The inclination to the ecliptic (i): "<<i_MA<<endl;
	cout<<"The longitude of the ascending node (Omega): "<<N_MA<<endl;
	cout<<"The argument of perihelion (w): "<<w_MA<<endl;
	cout<<"The mean anomaly (M): "<<M_MA<<endl;
	cout<<"The eccentric anomaly (E): "<<E_MA<<endl;
	cout<<"The true longitude (lambda): "<<lonecl_MA<<endl;
	cout<<"The true latitude (beta): "<<latecl_MA<<endl;
	cout<<"The distance (r): "<<r_MA<<endl;
	cout<<"Geocentric coordinates:"<<endl<<"x = "<<xg_MA<<endl<<"y = "<<yg_MA<<endl<<"z = "<<zg_MA<<endl;
	cout<<"Equatorial coordinates:"<<endl<<"x = "<<xe_MA<<endl<<"y = "<<ye_MA<<endl<<"z = "<<ze_MA<<endl<<"alpha = "<<RA_MA<<endl<<"delta = "<<Dec_MA<<endl;
	cout<<"Topocentric coordinates"<<endl<<"alpha = "<<topRA_MA<<endl<<"delta = "<<topDecl_MA<<endl;
	cout<<"*************************************************************************";

	cout<<endl<<"For Jupiter:"<<endl<<"*************************************************************************"<<endl;
	cout<<"The eccentricity (e): "<<e_J<<endl;
    cout<<"The mean distance from the earth (a): "<<a_J<<endl;
	cout<<"The inclination to the ecliptic (i): "<<i_J<<endl;
	cout<<"The longitude of the ascending node (Omega): "<<N_J<<endl;
	cout<<"The argument of perihelion (w): "<<w_J<<endl;
	cout<<"The mean anomaly (M): "<<M_J<<endl;
	cout<<"The eccentric anomaly (E): "<<E_J<<endl;
	cout<<"The true longitude (lambda): "<<lonecl_J<<endl;
	cout<<"The true latitude (beta): "<<latecl_J<<endl;
	cout<<"The distance (r): "<<r_J<<endl;
	cout<<"Geocentric coordinates:"<<endl<<"x = "<<xg_J<<endl<<"y = "<<yg_J<<endl<<"z = "<<zg_J<<endl;
	cout<<"Equatorial coordinates:"<<endl<<"x = "<<xe_J<<endl<<"y = "<<ye_J<<endl<<"z = "<<ze_J<<endl<<"alpha = "<<RA_J<<endl<<"delta = "<<Dec_J<<endl;
	cout<<"Topocentric coordinates"<<endl<<"alpha = "<<topRA_J<<endl<<"delta = "<<topDecl_J<<endl;
	cout<<"*************************************************************************";

	cout<<endl<<"For Saturn:"<<endl<<"*************************************************************************"<<endl;
	cout<<"The eccentricity (e): "<<e_SA<<endl;
    cout<<"The mean distance from the earth (a): "<<a_SA<<endl;
	cout<<"The inclination to the ecliptic (i): "<<i_SA<<endl;
	cout<<"The longitude of the ascending node (Omega): "<<N_SA<<endl;
	cout<<"The argument of perihelion (w): "<<w_SA<<endl;
	cout<<"The mean anomaly (M): "<<M_SA<<endl;
	cout<<"The eccentric anomaly (E): "<<E_SA<<endl;
	cout<<"The true longitude (lambda): "<<lonecl_SA<<endl;
	cout<<"The true latitude (beta): "<<latecl_SA<<endl;
	cout<<"The distance (r): "<<r_SA<<endl;
	cout<<"Geocentric coordinates:"<<endl<<"x = "<<xg_SA<<endl<<"y = "<<yg_SA<<endl<<"z = "<<zg_SA<<endl;
	cout<<"Equatorial coordinates:"<<endl<<"x = "<<xe_SA<<endl<<"y = "<<ye_SA<<endl<<"z = "<<ze_SA<<endl<<"alpha = "<<RA_SA<<endl<<"delta = "<<Dec_SA<<endl;
	cout<<"Topocentric coordinates"<<endl<<"alpha = "<<topRA_SA<<endl<<"delta = "<<topDecl_SA<<endl;
	cout<<"*************************************************************************";

	cout<<endl<<"For Uranus:"<<endl<<"*************************************************************************"<<endl;
	cout<<"The eccentricity (e): "<<e_U<<endl;
    cout<<"The mean distance from the earth (a): "<<a_U<<endl;
	cout<<"The inclination to the ecliptic (i): "<<i_U<<endl;
	cout<<"The longitude of the ascending node (Omega): "<<N_U<<endl;
	cout<<"The argument of perihelion (w): "<<w_U<<endl;
	cout<<"The mean anomaly (M): "<<M_U<<endl;
	cout<<"The eccentric anomaly (E): "<<E_U<<endl;
	cout<<"The true longitude (lambda): "<<lonecl_U<<endl;
	cout<<"The true latitude (beta): "<<latecl_U<<endl;
	cout<<"The distance (r): "<<r_U<<endl;
	cout<<"Geocentric coordinates:"<<endl<<"x = "<<xg_U<<endl<<"y = "<<yg_U<<endl<<"z = "<<zg_U<<endl;
	cout<<"Equatorial coordinates:"<<endl<<"x = "<<xe_U<<endl<<"y = "<<ye_U<<endl<<"z = "<<ze_U<<endl<<"alpha = "<<RA_U<<endl<<"delta = "<<Dec_U<<endl;
	cout<<"Topocentric coordinates"<<endl<<"alpha = "<<topRA_U<<endl<<"delta = "<<topDecl_U<<endl;
	cout<<"*************************************************************************";

	cout<<endl<<"For Neptune:"<<endl<<"*************************************************************************"<<endl;
	cout<<"The eccentricity (e): "<<e_N<<endl;
    cout<<"The mean distance from the earth (a): "<<a_N<<endl;
	cout<<"The inclination to the ecliptic (i): "<<i_N<<endl;
	cout<<"The longitude of the ascending node (Omega): "<<N_N<<endl;
	cout<<"The argument of perihelion (w): "<<w_N<<endl;
	cout<<"The mean anomaly (M): "<<M_N<<endl;
	cout<<"The eccentric anomaly (E): "<<E_N<<endl;
	cout<<"The true longitude (lambda): "<<lonecl_N<<endl;
	cout<<"The true latitude (beta): "<<latecl_N<<endl;
	cout<<"The distance (r): "<<r_N<<endl;
	cout<<"Geocentric coordinates:"<<endl<<"x = "<<xg_N<<endl<<"y = "<<yg_N<<endl<<"z = "<<zg_N<<endl;
	cout<<"Equatorial coordinates:"<<endl<<"x = "<<xe_N<<endl<<"y = "<<ye_N<<endl<<"z = "<<ze_N<<endl<<"alpha = "<<RA_N<<endl<<"delta = "<<Dec_N<<endl;
	cout<<"Topocentric coordinates"<<endl<<"alpha = "<<topRA_N<<endl<<"delta = "<<topDecl_N<<endl;
	cout<<"*************************************************************************";

	cout<<endl<<"For Pluto:"<<endl<<"*************************************************************************"<<endl;
	cout<<"The true longitude (lambda): "<<lonecl_P<<endl;
	cout<<"The true latitude (beta): "<<latecl_P<<endl;
	cout<<"The distance (r): "<<r_P<<endl;
	cout<<"Geocentric coordinates:"<<endl<<"x = "<<xg_P<<endl<<"y = "<<yg_P<<endl<<"z = "<<zg_P<<endl;
	cout<<"Equatorial coordinates:"<<endl<<"x = "<<xe_P<<endl<<"y = "<<ye_P<<endl<<"z = "<<ze_P<<endl<<"alpha = "<<RA_P<<endl<<"delta = "<<Dec_P<<endl;
	cout<<"Topocentric coordinates"<<endl<<"alpha = "<<topRA_P<<endl<<"delta = "<<topDecl_P<<endl;
	cout<<"*************************************************************************"<<endl<<endl;
	
	cout<<"Please enter the year:\n";
	cin>>y;
	cout<<"Please enter the month:\n";
	cin>>m;
	cout<<"Please enter the day:\n";
	cin>>Day;

	////////////////////////////////Set, rise and Meridian time////////////////////////////////
	d = 367*y - 7 * ( y + (m+9)/12) /4 + 275*m/9 + Day - 730530 ; //day
	const double reald = d;
	d = d + 0.3125; //LOCAL NOON
	a_S = 1.000000;
	ecl = 23.4393 - 3.563E-7 * d;  //obliquity of the ecliptic, i.e. the "tilt" of the Earth's axis of rotation
	e_S=0.016709 - 1.151E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    w_S=282.9404 + 4.70935E-5 * d;  //argument of perihelion
    M_S=356.0470 + 0.9856002585 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
	while(M_S>(360))
		M_S = M_S - 360;
	while(M_S<0)
		M_S = M_S + 360;
	L_S= M_S + w_S; //Mean Longitude of the Sun  (Ns=0)
    E_S = M_S + e_S*(180/pi) * sin(M_S*pi/180) * ( 1.0 + e_S * cos(M_S*pi/180) );  //Sun's eccentric anomaly E_S
    xv_S = a_S * (cos(E_S*pi/180) - e_S); //NOT IMPORTANT
    yv_S = a_S * sqrt(1.0 - e_S*e_S) * sin(E_S*pi/180); //NOT IMPORTANT
    vs = atan2( yv_S, xv_S ); //Sun's true anomaly
    rs = sqrt( xv_S*xv_S + yv_S*yv_S ); //Sun's distance
    lonsun = vs*180/pi + w_S; //sun's true longitude
	while(lonsun>360)
		lonsun = lonsun - 360;
    xs = rs * cos(lonsun*pi/180); //geocentric coordinates
    ys = rs * sin(lonsun*pi/180); //geocentric coordinates zs=0
    xes = xs; //equatorial coordinates
    yes = ys * cos(ecl*pi/180); //equatorial coordinates
    zes = ys * sin(ecl*pi/180); //equatorial coordinates
    RA_S  = atan2( yes, xes )*180/pi; //Sun's Right Ascension AT NOON LOCAL TIME
    Dec_S = atan2( zes, sqrt(xes*xes+yes*yes) )*180/pi; //Sun's Declination
	GMST0 = L_S + 180; //GSMT0 OF LOCAL NOON
	while(GMST0>(360))
		GMST0 = GMST0 - 360;
	while(GMST0<0)
		GMST0 = GMST0 + 360;
	double UT_Sun_in_south = ( RA_S - GMST0 - 51.5 );
	double cosLHA = (sin(-0.833*pi/180) - sin(35.0*pi/180) * sin(Dec_S*pi/180)) / (cos(35.0*pi/180) * cos(Dec_S*pi/180));
	double LHA;
	if(cosLHA<-1.0){
		cin>>D;
	}
	if(cosLHA>1.0){
		cin>>D;
	}
	double sunset;
	double sunrise;
	if(cosLHA<1 && cosLHA>-1)
	{
		LHA = acos(cosLHA)*180/pi; //Time in degrees
		if(LHA>(180))
			LHA = LHA - 360;
		if(LHA<-180)
			LHA = LHA + 360;
		sunset = UT_Sun_in_south + LHA;
		sunrise = UT_Sun_in_south - LHA;
	}
	while(sunset>(360))
		sunset = sunset - 360;
	while(sunset<0)
		sunset = sunset + 360;
	while(sunrise>(360))
		sunrise = sunrise - 360;
	while(sunrise<0)
		sunrise = sunrise + 360;
	sunset = sunset/15 + 4.5;
	sunrise = sunrise/15 + 4.5;
	while(UT_Sun_in_south>(360))
		UT_Sun_in_south = UT_Sun_in_south - 360;
	while(UT_Sun_in_south<0)
		UT_Sun_in_south = UT_Sun_in_south + 360;
	UT_Sun_in_south = UT_Sun_in_south/15 +4.5;
	double meridianhour = (int)UT_Sun_in_south;
	double meridianminute = (int)((UT_Sun_in_south - meridianhour) * 60);
	double meridiansecond = (int)((((UT_Sun_in_south - meridianhour) * 60) - meridianminute) * 60);
	double sunsethour = (int)sunset;
	double sunsetminute = (int)((sunset - sunsethour) * 60);
	double sunsetsecond = (int)((((sunset - sunsethour) * 60) - sunsetminute) * 60);
	double sunrisehour = (int)sunrise;
	double sunriseminute = (int)((sunrise - sunrisehour) * 60);
	double sunrisesecond = (int)((((sunrise - sunrisehour) * 60) - sunriseminute) * 60);
	cout<<"Sun:"<<endl<< "setting time: "<<sunsethour<<" : "<<sunsetminute<<" : "<<sunsetsecond<<endl<<"rising time: "<<sunrisehour<<" : "<<sunriseminute<<" : "<<sunrisesecond<<endl<<"meridian time:"<<meridianhour<<" : "<<meridianminute<<" : "<<meridiansecond<<endl;;
	
	////////////////////////////////Moon's rise and set times////////////////////////////////
	 N_M=125.1228 - 0.0529538083 * d;  //longitude of the ascending node
	while(N_M>360)
		N_M = N_M - 360;
	while(N_M<0)
		N_M = N_M + 360;
     i_M=5.1454;  //inclination to the ecliptic (plane of the Earth's orbit)
     w_M=318.0634 + 0.1643573223 * d;  //argument of perihelion
	while(w_M>(360))
		w_M = w_M - 360;
	while(w_M<0)
		w_M = w_M + 360;
     a_M=60.2666;  //semi-major axis, or mean distance from Sun (earth radii)
     e_M=0.054900;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
     M_M=115.3654 + 13.0649929509 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_M>(360))
		M_M = M_M - 360;
	while(M_M<0)
		M_M = M_M + 360;
	 L_M = M_M + w_M + N_M ; //Mean longitude of the Moon
	 D= L_M - L_S; //Mean elongation of the Moon
	 F= L_M - N_M; //Argument of latitude for the Moon
    
	  E_M = M_M + e_M*(180/pi) * sin(M_M*pi/180) * ( 1.0 + e_M * cos(M_M*pi/180) );  //Moon's eccentric anomaly E_S
    if(e_M>0.05){
        double E0_M=E_M;
        double E1_M=E_M + 2;
        while(fabs(E0_M-E1_M)>0.001){
            E0_M=E1_M;
            E1_M = E_M - ( E_M - e_M*(180/pi) * sin(E_M*pi/180) - M_M ) / ( 1 - e_M * cos(E_M*pi/180) );
        }
		E_M=E1_M;
    }
     xv_M = a_M * (cos(E_M*pi/180) - e_M); //NOT IMPORTANT
     yv_M = a_M * sqrt(1.0 - e_M*e_M) * sin(E_M*pi/180); //NOT IMPORTANT
     v_M = atan2( yv_M, xv_M ); //Moon's true anomaly
     r_M = sqrt( xv_M*xv_M + yv_M*yv_M ); //Moon's distance
     xh_M = r_M * ( cos(N_M*pi/180) * cos(v_M+w_M*pi/180) - sin(N_M*pi/180) * sin(v_M+w_M*pi/180) * cos(i_M*pi/180) ); //Moon's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_M = r_M * ( sin(N_M*pi/180) * cos(v_M+w_M*pi/180) + cos(N_M*pi/180) * sin(v_M+w_M*pi/180) * cos(i_M*pi/180) ); //Moon's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_M = r_M * ( sin(v_M+w_M*pi/180) * sin(i_M*pi/180) ); //Moon's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_M = atan2( yh_M, xh_M )*180/pi; //Moon's ecliptic longitude
     latecl_M = atan2( zh_M, sqrt(xh_M*xh_M+yh_M*yh_M) )*180/pi; //Moon's ecliptic latitude
	
	lonecl_M = lonecl_M -1.274 * sin(M_M*pi/180 - 2*D*pi/180)          
    +0.658 * sin(2*D*pi/180)               
    -0.186 * sin(M_S*pi/180)                
    -0.059 * sin(2*M_M*pi/180 - 2*D*pi/180)
    -0.057 * sin(M_M*pi/180 - 2*D*pi/180 + M_S*pi/180)
    +0.053 * sin(M_M*pi/180 + 2*D*pi/180)
    +0.046 * sin(2*D*pi/180 - M_S*pi/180)
    +0.041 * sin(M_M*pi/180 - M_S*pi/180)
    -0.035 * sin(D*pi/180)                 
    -0.031 * sin(M_M*pi/180 + M_S*pi/180)
    -0.015 * sin(2*F*pi/180 - 2*D*pi/180)
    +0.011 * sin(M_M*pi/180 - 4*D*pi/180);

	latecl_M = latecl_M -0.173 * sin(F*pi/180 - 2*D*pi/180)
    -0.055 * sin(M_M*pi/180 - F*pi/180 - 2*D*pi/180)
    -0.046 * sin(M_M*pi/180 + F*pi/180 - 2*D*pi/180)
    +0.033 * sin(F*pi/180 + 2*D*pi/180)
    +0.017 * sin(2*M_M*pi/180 + F*pi/180);

	r_M = r_M -0.58 * cos(M_M*pi/180 - 2*D*pi/180)
    -0.46 * cos(2*D*pi/180);

	 xg_M = r_M * cos(lonecl_M*pi/180) * cos(latecl_M*pi/180); //xh_M + xs; //geocentric position
	 yg_M = r_M * sin(lonecl_M*pi/180) * cos(latecl_M*pi/180); //yh_M + ys; //geocentric position
	 zg_M = r_M * sin(latecl_M*pi/180); //zh_M; //geocentric position
	 xe_M = xg_M; //equatorial coordinates
     ye_M = yg_M * cos(ecl*pi/180) - zg_M * sin(ecl*pi/180); //equatorial coordinates
     ze_M = yg_M * sin(ecl*pi/180) + zg_M * cos(ecl*pi/180); //equatorial coordinates
     RA_M  = atan2( ye_M, xe_M )*180/pi; //Right Ascension
     Dec_M = atan2( ze_M, sqrt(xe_M*xe_M+ye_M*ye_M) )*180/pi; //Declination
	 HA_M = LST - RA_M; //hour angle of the Moon
	 x_M = cos(HA_M*pi/180) * cos(Dec_M*pi/180);
     y_M = sin(HA_M*pi/180) * cos(Dec_M*pi/180);
     z_M = sin(Dec_M*pi/180);
	 xhor_M = x_M * sin(35*pi/180) - z_M * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_M = y_M; //Horizontal coordinates
     zhor_M = x_M * cos(35*pi/180) + z_M * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_M  = atan2( yhor_M, xhor_M )*180/pi + 180; //local azimuth
     alt_M = atan2( zhor_M, sqrt(xhor_M*xhor_M+yhor_M*yhor_M) )*180/pi; //local altitude
	 mpar = asin( 1/r_M )*180/pi; //Moon's Paralax
	 alt_topoc_M = alt_M - mpar * cos(alt_M*pi/180); //topocentric altitude
	 gclat_M = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_M   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	//double gclat = 35;
	//double rho = 1.0;
	 g_M = atan( tan(gclat_M*pi/180) / cos(HA_M*pi/180) ); //an auxiliary angle
	 topRA_M   = RA_M  - mpar * rho_M * cos(gclat_M*pi/180) * sin(HA_M*pi/180) / cos(Dec_M*pi/180);
     topDecl_M = Dec_M - mpar * rho_M * sin(gclat_M*pi/180) * sin(g_M - Dec_M*pi/180) / sin(g_M);

	 double UT_M_in_south = ( topRA_M - GMST0 - 51.5 );
	while(UT_M_in_south>360)
		UT_M_in_south = UT_M_in_south - 360;
	while(UT_M_in_south<0)
		UT_M_in_south = UT_M_in_south +360;
	double cosLHA_M = (sin((-0.583-0.266667)*pi/180) - sin(35.0*pi/180) * sin(topDecl_M*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_M*pi/180));
	double LHA_M;
	if(cosLHA_M<-1.0){
		cin>>D;
	}
	if(cosLHA_M>1.0){
		cin>>D;
	}
	double set_M;
	double rise_M;
	if(cosLHA_M<1 && cosLHA_M>-1)
	{
		LHA_M = acos(cosLHA_M)*180/pi; //Time in degrees
		if(LHA_M>(180))
			LHA_M = LHA_M- 360;
		if(LHA_M<-180)
			LHA_M = LHA_M + 360;
		set_M = UT_M_in_south + LHA_M;
		rise_M = UT_M_in_south - LHA_M;
	}
	while(set_M>(360))
		set_M = set_M - 360;
	while(set_M<0)
		set_M = set_M + 360;
	while(rise_M>(360))
		rise_M = rise_M - 360;
	while(rise_M<0)
		rise_M = rise_M + 360;
	set_M = set_M/15.04107  +4.5;
	rise_M = rise_M/15.04107  +4.5;
	//////////////////////////////////////////////////////////////////////////////////PERTURBING SET
	double perturbed_rise_time_M = rise_M/24;
	double perturbed_set_time_M = set_M/24;


	for(int set=0;set<4;set++){
	d = 367*y - 7 * ( y + (m+9)/12) /4 + 275*m/9 + Day - 730530 ; //day
	d = d + perturbed_set_time_M; //perturbing time
		 N_M=125.1228 - 0.0529538083 * d;  //longitude of the ascending node
	while(N_M>360)
		N_M = N_M - 360;
	while(N_M<0)
		N_M = N_M + 360;
     i_M=5.1454;  //inclination to the ecliptic (plane of the Earth's orbit)
     w_M=318.0634 + 0.1643573223 * d;  //argument of perihelion
	while(w_M>(360))
		w_M = w_M - 360;
	while(w_M<0)
		w_M = w_M + 360;
     a_M=60.2666;  //semi-major axis, or mean distance from Sun (earth radii)
     e_M=0.054900;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
     M_M=115.3654 + 13.0649929509 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_M>(360))
		M_M = M_M - 360;
	while(M_M<0)
		M_M = M_M + 360;
	 L_M = M_M + w_M + N_M ; //Mean longitude of the Moon
	 D= L_M - L_S; //Mean elongation of the Moon
	 F= L_M - N_M; //Argument of latitude for the Moon
    
	  E_M = M_M + e_M*(180/pi) * sin(M_M*pi/180) * ( 1.0 + e_M * cos(M_M*pi/180) );  //Moon's eccentric anomaly E_S
    if(e_M>0.05){
        double E0_M=E_M;
        double E1_M=E_M + 2;
        while(fabs(E0_M-E1_M)>0.001){
            E0_M=E1_M;
            E1_M = E_M - ( E_M - e_M*(180/pi) * sin(E_M*pi/180) - M_M ) / ( 1 - e_M * cos(E_M*pi/180) );
        }
		E_M=E1_M;
    }
     xv_M = a_M * (cos(E_M*pi/180) - e_M); //NOT IMPORTANT
     yv_M = a_M * sqrt(1.0 - e_M*e_M) * sin(E_M*pi/180); //NOT IMPORTANT
     v_M = atan2( yv_M, xv_M ); //Moon's true anomaly
     r_M = sqrt( xv_M*xv_M + yv_M*yv_M ); //Moon's distance
     xh_M = r_M * ( cos(N_M*pi/180) * cos(v_M+w_M*pi/180) - sin(N_M*pi/180) * sin(v_M+w_M*pi/180) * cos(i_M*pi/180) ); //Moon's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_M = r_M * ( sin(N_M*pi/180) * cos(v_M+w_M*pi/180) + cos(N_M*pi/180) * sin(v_M+w_M*pi/180) * cos(i_M*pi/180) ); //Moon's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_M = r_M * ( sin(v_M+w_M*pi/180) * sin(i_M*pi/180) ); //Moon's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_M = atan2( yh_M, xh_M )*180/pi; //Moon's ecliptic longitude
     latecl_M = atan2( zh_M, sqrt(xh_M*xh_M+yh_M*yh_M) )*180/pi; //Moon's ecliptic latitude
	
	lonecl_M = lonecl_M -1.274 * sin(M_M*pi/180 - 2*D*pi/180)          
    +0.658 * sin(2*D*pi/180)               
    -0.186 * sin(M_S*pi/180)                
    -0.059 * sin(2*M_M*pi/180 - 2*D*pi/180)
    -0.057 * sin(M_M*pi/180 - 2*D*pi/180 + M_S*pi/180)
    +0.053 * sin(M_M*pi/180 + 2*D*pi/180)
    +0.046 * sin(2*D*pi/180 - M_S*pi/180)
    +0.041 * sin(M_M*pi/180 - M_S*pi/180)
    -0.035 * sin(D*pi/180)                 
    -0.031 * sin(M_M*pi/180 + M_S*pi/180)
    -0.015 * sin(2*F*pi/180 - 2*D*pi/180)
    +0.011 * sin(M_M*pi/180 - 4*D*pi/180);

	latecl_M = latecl_M -0.173 * sin(F*pi/180 - 2*D*pi/180)
    -0.055 * sin(M_M*pi/180 - F*pi/180 - 2*D*pi/180)
    -0.046 * sin(M_M*pi/180 + F*pi/180 - 2*D*pi/180)
    +0.033 * sin(F*pi/180 + 2*D*pi/180)
    +0.017 * sin(2*M_M*pi/180 + F*pi/180);

	r_M = r_M -0.58 * cos(M_M*pi/180 - 2*D*pi/180)
    -0.46 * cos(2*D*pi/180);

	 xg_M = r_M * cos(lonecl_M*pi/180) * cos(latecl_M*pi/180); //xh_M + xs; //geocentric position
	 yg_M = r_M * sin(lonecl_M*pi/180) * cos(latecl_M*pi/180); //yh_M + ys; //geocentric position
	 zg_M = r_M * sin(latecl_M*pi/180); //zh_M; //geocentric position
	 xe_M = xg_M; //equatorial coordinates
     ye_M = yg_M * cos(ecl*pi/180) - zg_M * sin(ecl*pi/180); //equatorial coordinates
     ze_M = yg_M * sin(ecl*pi/180) + zg_M * cos(ecl*pi/180); //equatorial coordinates
     RA_M  = atan2( ye_M, xe_M )*180/pi; //Right Ascension
     Dec_M = atan2( ze_M, sqrt(xe_M*xe_M+ye_M*ye_M) )*180/pi; //Declination
	 HA_M = LST - RA_M; //hour angle of the Moon
	 x_M = cos(HA_M*pi/180) * cos(Dec_M*pi/180);
     y_M = sin(HA_M*pi/180) * cos(Dec_M*pi/180);
     z_M = sin(Dec_M*pi/180);
	 xhor_M = x_M * sin(35*pi/180) - z_M * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_M = y_M; //Horizontal coordinates
     zhor_M = x_M * cos(35*pi/180) + z_M * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_M  = atan2( yhor_M, xhor_M )*180/pi + 180; //local azimuth
     alt_M = atan2( zhor_M, sqrt(xhor_M*xhor_M+yhor_M*yhor_M) )*180/pi; //local altitude
	 mpar = asin( 1/r_M )*180/pi; //Moon's Paralax
	 alt_topoc_M = alt_M - mpar * cos(alt_M*pi/180); //topocentric altitude
	 gclat_M = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_M   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	//double gclat = 35;
	//double rho = 1.0;
	 g_M = atan( tan(gclat_M*pi/180) / cos(HA_M*pi/180) ); //an auxiliary angle
	 topRA_M   = RA_M  - mpar * rho_M * cos(gclat_M*pi/180) * sin(HA_M*pi/180) / cos(Dec_M*pi/180);
     topDecl_M = Dec_M - mpar * rho_M * sin(gclat_M*pi/180) * sin(g_M - Dec_M*pi/180) / sin(g_M);

	  UT_M_in_south = ( topRA_M - GMST0 - 51.5 );
	while(UT_M_in_south>360)
		UT_M_in_south = UT_M_in_south - 360;
	while(UT_M_in_south<0)
		UT_M_in_south = UT_M_in_south +360;
	 cosLHA_M = (sin((-0.583)*pi/180) - sin(35.0*pi/180) * sin(topDecl_M*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_M*pi/180));

	if(cosLHA_M<-1.0){
		cin>>D;
	}
	if(cosLHA_M>1.0){
		cin>>D;
	}
	 set_M;
	rise_M;
	if(cosLHA_M<1 && cosLHA_M>-1)
	{
		LHA_M = acos(cosLHA_M)*180/pi; //Time in degrees
		if(LHA_M>(180))
			LHA_M = LHA_M- 360;
		if(LHA_M<-180)
			LHA_M = LHA_M + 360;
		set_M = UT_M_in_south + LHA_M;
	}
	while(set_M>(360))
		set_M = set_M - 360;
	while(set_M<0)
		set_M = set_M + 360;
	set_M = set_M/15.04107  +4.5;
	perturbed_set_time_M = set_M/24;
	}
	perturbed_set_time_M = perturbed_set_time_M*24;

	
	//////////////////////////////////////////////////////////////////////////////////PERTURBING RISE

	for(int rise=0;rise<4;rise++){
	d = 367*y - 7 * ( y + (m+9)/12) /4 + 275*m/9 + Day - 730530 ; //day
	d = d + perturbed_rise_time_M; //perturbing time
		 N_M=125.1228 - 0.0529538083 * d;  //longitude of the ascending node
	while(N_M>360)
		N_M = N_M - 360;
	while(N_M<0)
		N_M = N_M + 360;
     i_M=5.1454;  //inclination to the ecliptic (plane of the Earth's orbit)
     w_M=318.0634 + 0.1643573223 * d;  //argument of perihelion
	while(w_M>(360))
		w_M = w_M - 360;
	while(w_M<0)
		w_M = w_M + 360;
     a_M=60.2666;  //semi-major axis, or mean distance from Sun (earth radii)
     e_M=0.054900;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
     M_M=115.3654 + 13.0649929509 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_M>(360))
		M_M = M_M - 360;
	while(M_M<0)
		M_M = M_M + 360;
	 L_M = M_M + w_M + N_M ; //Mean longitude of the Moon
	 D= L_M - L_S; //Mean elongation of the Moon
	 F= L_M - N_M; //Argument of latitude for the Moon
    
	  E_M = M_M + e_M*(180/pi) * sin(M_M*pi/180) * ( 1.0 + e_M * cos(M_M*pi/180) );  //Moon's eccentric anomaly E_S
    if(e_M>0.05){
        double E0_M=E_M;
        double E1_M=E_M + 2;
        while(fabs(E0_M-E1_M)>0.001){
            E0_M=E1_M;
            E1_M = E_M - ( E_M - e_M*(180/pi) * sin(E_M*pi/180) - M_M ) / ( 1 - e_M * cos(E_M*pi/180) );
        }
		E_M=E1_M;
    }
     xv_M = a_M * (cos(E_M*pi/180) - e_M); //NOT IMPORTANT
     yv_M = a_M * sqrt(1.0 - e_M*e_M) * sin(E_M*pi/180); //NOT IMPORTANT
     v_M = atan2( yv_M, xv_M ); //Moon's true anomaly
     r_M = sqrt( xv_M*xv_M + yv_M*yv_M ); //Moon's distance
     xh_M = r_M * ( cos(N_M*pi/180) * cos(v_M+w_M*pi/180) - sin(N_M*pi/180) * sin(v_M+w_M*pi/180) * cos(i_M*pi/180) ); //Moon's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_M = r_M * ( sin(N_M*pi/180) * cos(v_M+w_M*pi/180) + cos(N_M*pi/180) * sin(v_M+w_M*pi/180) * cos(i_M*pi/180) ); //Moon's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_M = r_M * ( sin(v_M+w_M*pi/180) * sin(i_M*pi/180) ); //Moon's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_M = atan2( yh_M, xh_M )*180/pi; //Moon's ecliptic longitude
     latecl_M = atan2( zh_M, sqrt(xh_M*xh_M+yh_M*yh_M) )*180/pi; //Moon's ecliptic latitude
	
	lonecl_M = lonecl_M -1.274 * sin(M_M*pi/180 - 2*D*pi/180)          
    +0.658 * sin(2*D*pi/180)               
    -0.186 * sin(M_S*pi/180)                
    -0.059 * sin(2*M_M*pi/180 - 2*D*pi/180)
    -0.057 * sin(M_M*pi/180 - 2*D*pi/180 + M_S*pi/180)
    +0.053 * sin(M_M*pi/180 + 2*D*pi/180)
    +0.046 * sin(2*D*pi/180 - M_S*pi/180)
    +0.041 * sin(M_M*pi/180 - M_S*pi/180)
    -0.035 * sin(D*pi/180)                 
    -0.031 * sin(M_M*pi/180 + M_S*pi/180)
    -0.015 * sin(2*F*pi/180 - 2*D*pi/180)
    +0.011 * sin(M_M*pi/180 - 4*D*pi/180);

	latecl_M = latecl_M -0.173 * sin(F*pi/180 - 2*D*pi/180)
    -0.055 * sin(M_M*pi/180 - F*pi/180 - 2*D*pi/180)
    -0.046 * sin(M_M*pi/180 + F*pi/180 - 2*D*pi/180)
    +0.033 * sin(F*pi/180 + 2*D*pi/180)
    +0.017 * sin(2*M_M*pi/180 + F*pi/180);

	r_M = r_M -0.58 * cos(M_M*pi/180 - 2*D*pi/180)
    -0.46 * cos(2*D*pi/180);

	 xg_M = r_M * cos(lonecl_M*pi/180) * cos(latecl_M*pi/180); //xh_M + xs; //geocentric position
	 yg_M = r_M * sin(lonecl_M*pi/180) * cos(latecl_M*pi/180); //yh_M + ys; //geocentric position
	 zg_M = r_M * sin(latecl_M*pi/180); //zh_M; //geocentric position
	 xe_M = xg_M; //equatorial coordinates
     ye_M = yg_M * cos(ecl*pi/180) - zg_M * sin(ecl*pi/180); //equatorial coordinates
     ze_M = yg_M * sin(ecl*pi/180) + zg_M * cos(ecl*pi/180); //equatorial coordinates
     RA_M  = atan2( ye_M, xe_M )*180/pi; //Right Ascension
     Dec_M = atan2( ze_M, sqrt(xe_M*xe_M+ye_M*ye_M) )*180/pi; //Declination
	 HA_M = LST - RA_M; //hour angle of the Moon
	 x_M = cos(HA_M*pi/180) * cos(Dec_M*pi/180);
     y_M = sin(HA_M*pi/180) * cos(Dec_M*pi/180);
     z_M = sin(Dec_M*pi/180);
	 xhor_M = x_M * sin(35*pi/180) - z_M * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_M = y_M; //Horizontal coordinates
     zhor_M = x_M * cos(35*pi/180) + z_M * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_M  = atan2( yhor_M, xhor_M )*180/pi + 180; //local azimuth
     alt_M = atan2( zhor_M, sqrt(xhor_M*xhor_M+yhor_M*yhor_M) )*180/pi; //local altitude
	 mpar = asin( 1/r_M )*180/pi; //Moon's Paralax
	 alt_topoc_M = alt_M - mpar * cos(alt_M*pi/180); //topocentric altitude
	 gclat_M = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_M   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	//double gclat = 35;
	//double rho = 1.0;
	 g_M = atan( tan(gclat_M*pi/180) / cos(HA_M*pi/180) ); //an auxiliary angle
	 topRA_M   = RA_M  - mpar * rho_M * cos(gclat_M*pi/180) * sin(HA_M*pi/180) / cos(Dec_M*pi/180);
     topDecl_M = Dec_M - mpar * rho_M * sin(gclat_M*pi/180) * sin(g_M - Dec_M*pi/180) / sin(g_M);

	  UT_M_in_south = ( topRA_M - GMST0 - 51.5 );
	while(UT_M_in_south>360)
		UT_M_in_south = UT_M_in_south - 360;
	while(UT_M_in_south<0)
		UT_M_in_south = UT_M_in_south +360;
	 cosLHA_M = (sin(-0.583*pi/180) - sin(35.0*pi/180) * sin(topDecl_M*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_M*pi/180));

	if(cosLHA_M<-1.0){
		cin>>D;
	}
	if(cosLHA_M>1.0){
		cin>>D;
	}
	if(cosLHA_M<1 && cosLHA_M>-1)
	{
		LHA_M = acos(cosLHA_M)*180/pi; //Time in degrees
		if(LHA_M>(180))
			LHA_M = LHA_M- 360;
		if(LHA_M<-180)
			LHA_M = LHA_M + 360;
		rise_M = UT_M_in_south - LHA_M;
	}
	while(rise_M>(360))
		rise_M = rise_M - 360;
	while(rise_M<0)
		rise_M = rise_M + 360;
	rise_M = rise_M/15.04107  + 4.5;
	perturbed_rise_time_M = rise_M/24;
	}
	perturbed_rise_time_M = perturbed_rise_time_M*24;
	UT_M_in_south = UT_M_in_south/15;
	UT_M_in_south = UT_M_in_south + 4.5;
	double meridianhour_M = (int)UT_M_in_south;
	double meridianminute_M = (int)((UT_M_in_south - meridianhour_M) * 60);
	double meridiansecond_M = (int)((((UT_M_in_south - meridianhour_M) * 60) - meridianminute_M) * 60);
	double sethour_M = (int)perturbed_set_time_M;
	double setminute_M = (int)((perturbed_set_time_M - sethour_M) * 60);
	double setsecond_M = (int)((((perturbed_set_time_M - sethour_M) * 60) - setminute_M) * 60);
	double risehour_M = (int)perturbed_rise_time_M;
	double riseminute_M = (int)((perturbed_rise_time_M - risehour_M) * 60);
	double risesecond_M = (int)((((perturbed_rise_time_M - risehour_M) * 60) - riseminute_M) * 60);
	cout<<"Moon:"<<endl<< "setting time: "<<sethour_M<<" : "<<setminute_M<<" : "<<setsecond_M<<endl<<"rising time: "<<risehour_M<<" : "<<riseminute_M<<" : "<<risesecond_M<<endl<<"meridian time:"<<meridianhour_M<<" : "<<meridianminute_M<<" : "<<meridiansecond_M<<endl;

	////////////////////////////////Mercury's rise and set times//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	d = reald + 0.3125;
	N_ME=48.3313 + 3.24587E-5 * d;  //longitude of the ascending node
    i_ME=7.0047 + 5.00E-8 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    w_ME=29.1241 + 1.01444E-5 * d;  //argument of perihelion
    a_ME=0.387098;  //semi-major axis, or mean distance from Sun
    e_ME=0.205635 + 5.59E-10 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M_ME=168.6562 + 4.0923344368 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_ME>(360))
		M_ME = M_ME - 360;
	while(M_ME<0)
		M_ME = M_ME + 360;
    
	  E_ME = M_ME + e_ME*(180/pi) * sin(M_ME*pi/180) * ( 1.0 + e_ME * cos(M_ME*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_ME>0.05){
        double E0_ME=E_ME;
        double E1_ME=E_ME;
        while(fabs(E0_ME-E1_ME)>0.001){
            E0_ME=E1_ME;
            E1_ME = E_ME - ( E_ME - e_ME*(180/pi) * sin(E_ME*pi/180) - M_ME ) / ( 1 - e_ME * cos(E_ME*pi/180) );
        }
    }
     xv_ME = a_ME * (cos(E_ME*pi/180) - e_ME); //NOT IMPORTANT
     yv_ME = a_ME * sqrt(1.0 - e_ME*e_ME) * sin(E_ME*pi/180); //NOT IMPORTANT
     v_ME = atan2( yv_ME, xv_ME ); //Mercury's true anomaly
     r_ME = sqrt( xv_ME*xv_ME + yv_ME*yv_ME ); //Mercury's distance
     xh_ME = r_ME * ( cos(N_ME*pi/180) * cos(v_ME+w_ME*pi/180) - sin(N_ME*pi/180) * sin(v_ME+w_ME*pi/180) * cos(i_ME*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_ME = r_ME * ( sin(N_ME*pi/180) * cos(v_ME+w_ME*pi/180) + cos(N_ME*pi/180) * sin(v_ME+w_ME*pi/180) * cos(i_ME*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_ME = r_ME * ( sin(v_ME+w_ME*pi/180) * sin(i_ME*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_ME = atan2( yh_ME, xh_ME )*180/pi; //Mercury's ecliptic longitude
     latecl_ME = atan2( zh_ME, sqrt(xh_ME*xh_ME+yh_ME*yh_ME) )*180/pi; //Mercury's ecliptic latitude
	 xg_ME = r_ME * cos(lonecl_ME*pi/180) * cos(latecl_ME*pi/180) +xs; //geocentric position
	 yg_ME = r_ME * sin(lonecl_ME*pi/180) * cos(latecl_ME*pi/180) +ys; //geocentric position
	 zg_ME = r_ME * sin(latecl_ME*pi/180); //geocentric position
	 xe_ME = xg_ME; //equatorial coordinates
     ye_ME = yg_ME * cos(ecl*pi/180) - zg_ME * sin(ecl*pi/180); //equatorial coordinates
     ze_ME = yg_ME * sin(ecl*pi/180) + zg_ME * cos(ecl*pi/180); //equatorial coordinates
     RA_ME  = atan2( ye_ME, xe_ME )*180/pi; //Right Ascension
     Dec_ME = atan2( ze_ME, sqrt(xe_ME*xe_ME+ye_ME*ye_ME) )*180/pi; //Declination
	 par_ME = (8.794/3600) / r_ME;
	 HA_ME = LST - RA_ME; //hour angle of the Moon
	 x_ME = cos(HA_ME*pi/180) * cos(Dec_ME*pi/180);
     y_ME = sin(HA_ME*pi/180) * cos(Dec_ME*pi/180);
     z_ME = sin(Dec_ME*pi/180);
	 xhor_ME = x_ME * sin(35*pi/180) - z_ME * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_ME = y_ME; //Horizontal coordinates
     zhor_ME = x_ME * cos(35*pi/180) + z_ME * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_ME  = atan2( yhor_ME, xhor_ME )*180/pi + 180; //local azimuth
     alt_ME = atan2( zhor_ME, sqrt(xhor_ME*xhor_ME+yhor_ME*yhor_ME) )*180/pi; //local altitude
	 alt_topoc_ME = alt_ME - par_ME * cos(alt_ME*pi/180); //topocentric altitude
	 gclat_ME = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_ME   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_ME = atan( tan(gclat_ME*pi/180) / cos(HA_ME*pi/180) ); //an auxiliary angle
	 topRA_ME   = RA_ME  - par_ME * rho_ME * cos(gclat_ME*pi/180) * sin(HA_ME*pi/180) / cos(Dec_ME*pi/180);
     topDecl_ME = Dec_ME - par_ME * rho_ME * sin(gclat_ME*pi/180) * sin(g_ME - Dec_ME*pi/180) / sin(g_ME);
	
	 double UT_ME_in_south = ( topRA_ME - GMST0 - 51.5 );
	while(UT_ME_in_south>360)
		UT_ME_in_south = UT_ME_in_south - 360;
	while(UT_ME_in_south<0)
		UT_ME_in_south = UT_ME_in_south +360;
	double cosLHA_ME = (sin((-0.583)*pi/180) - sin(35.0*pi/180) * sin(topDecl_ME*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_ME*pi/180));
	double LHA_ME;
	if(cosLHA_ME<-1.0){
		cin>>D;
	}
	if(cosLHA_ME>1.0){
		cin>>D;
	}
	double set_ME;
	double rise_ME;
	if(cosLHA_ME<1 && cosLHA_ME>-1)
	{
		LHA_ME = acos(cosLHA_ME)*180/pi; //Time in degrees
		if(LHA_ME>(180))
			LHA_ME = LHA_ME- 360;
		if(LHA_ME<-180)
			LHA_ME = LHA_ME + 360;
		set_ME = UT_ME_in_south + LHA_ME;
		rise_ME = UT_ME_in_south - LHA_ME;
	}
	while(set_ME>(360))
		set_ME = set_ME - 360;
	while(set_ME<0)
		set_ME = set_ME + 360;
	while(rise_ME>(360))
		rise_ME = rise_ME - 360;
	while(rise_ME<0)
		rise_ME = rise_ME + 360;
	set_ME = set_ME/15.04107  +4.5;
	rise_ME = rise_ME/15.04107  +4.5;
	//////////////////////////////////////////////////////////////////////////////////PERTURBING SET
	double perturbed_rise_time_ME = rise_ME/24;
	double perturbed_set_time_ME = set_ME/24;

	for(int set=0;set<4;set++){
	d = reald ; //day
	d = d + perturbed_set_time_ME; //perturbing time
	N_ME=48.3313 + 3.24587E-5 * d;  //longitude of the ascending node
    i_ME=7.0047 + 5.00E-8 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    w_ME=29.1241 + 1.01444E-5 * d;  //argument of perihelion
    a_ME=0.387098;  //semi-major axis, or mean distance from Sun
    e_ME=0.205635 + 5.59E-10 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M_ME=168.6562 + 4.0923344368 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_ME>(360))
		M_ME = M_ME - 360;
	while(M_ME<0)
		M_ME = M_ME + 360;
    
	  E_ME = M_ME + e_ME*(180/pi) * sin(M_ME*pi/180) * ( 1.0 + e_ME * cos(M_ME*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_ME>0.05){
        double E0_ME=E_ME;
        double E1_ME=E_ME;
        while(fabs(E0_ME-E1_ME)>0.001){
            E0_ME=E1_ME;
            E1_ME = E_ME - ( E_ME - e_ME*(180/pi) * sin(E_ME*pi/180) - M_ME ) / ( 1 - e_ME * cos(E_ME*pi/180) );
        }
    }
     xv_ME = a_ME * (cos(E_ME*pi/180) - e_ME); //NOT IMPORTANT
     yv_ME = a_ME * sqrt(1.0 - e_ME*e_ME) * sin(E_ME*pi/180); //NOT IMPORTANT
     v_ME = atan2( yv_ME, xv_ME ); //Mercury's true anomaly
     r_ME = sqrt( xv_ME*xv_ME + yv_ME*yv_ME ); //Mercury's distance
     xh_ME = r_ME * ( cos(N_ME*pi/180) * cos(v_ME+w_ME*pi/180) - sin(N_ME*pi/180) * sin(v_ME+w_ME*pi/180) * cos(i_ME*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_ME = r_ME * ( sin(N_ME*pi/180) * cos(v_ME+w_ME*pi/180) + cos(N_ME*pi/180) * sin(v_ME+w_ME*pi/180) * cos(i_ME*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_ME = r_ME * ( sin(v_ME+w_ME*pi/180) * sin(i_ME*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_ME = atan2( yh_ME, xh_ME )*180/pi; //Mercury's ecliptic longitude
     latecl_ME = atan2( zh_ME, sqrt(xh_ME*xh_ME+yh_ME*yh_ME) )*180/pi; //Mercury's ecliptic latitude
	 xg_ME = r_ME * cos(lonecl_ME*pi/180) * cos(latecl_ME*pi/180) +xs; //geocentric position
	 yg_ME = r_ME * sin(lonecl_ME*pi/180) * cos(latecl_ME*pi/180) +ys; //geocentric position
	 zg_ME = r_ME * sin(latecl_ME*pi/180); //geocentric position
	 xe_ME = xg_ME; //equatorial coordinates
     ye_ME = yg_ME * cos(ecl*pi/180) - zg_ME * sin(ecl*pi/180); //equatorial coordinates
     ze_ME = yg_ME * sin(ecl*pi/180) + zg_ME * cos(ecl*pi/180); //equatorial coordinates
     RA_ME  = atan2( ye_ME, xe_ME )*180/pi; //Right Ascension
     Dec_ME = atan2( ze_ME, sqrt(xe_ME*xe_ME+ye_ME*ye_ME) )*180/pi; //Declination
	 par_ME = (8.794/3600) / r_ME;
	 HA_ME = LST - RA_ME; //hour angle of the Moon
	 x_ME = cos(HA_ME*pi/180) * cos(Dec_ME*pi/180);
     y_ME = sin(HA_ME*pi/180) * cos(Dec_ME*pi/180);
     z_ME = sin(Dec_ME*pi/180);
	 xhor_ME = x_ME * sin(35*pi/180) - z_ME * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_ME = y_ME; //Horizontal coordinates
     zhor_ME = x_ME * cos(35*pi/180) + z_ME * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_ME  = atan2( yhor_ME, xhor_ME )*180/pi + 180; //local azimuth
     alt_ME = atan2( zhor_ME, sqrt(xhor_ME*xhor_ME+yhor_ME*yhor_ME) )*180/pi; //local altitude
	 alt_topoc_ME = alt_ME - par_ME * cos(alt_ME*pi/180); //topocentric altitude
	 gclat_ME = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_ME   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_ME = atan( tan(gclat_ME*pi/180) / cos(HA_ME*pi/180) ); //an auxiliary angle
	 topRA_ME   = RA_ME  - par_ME * rho_ME * cos(gclat_ME*pi/180) * sin(HA_ME*pi/180) / cos(Dec_ME*pi/180);
     topDecl_ME = Dec_ME - par_ME * rho_ME * sin(gclat_ME*pi/180) * sin(g_ME - Dec_ME*pi/180) / sin(g_ME);
	 
	  UT_ME_in_south = ( topRA_ME - GMST0 - 51.5 );
	while(UT_ME_in_south>360)
		UT_ME_in_south = UT_ME_in_south - 360;
	while(UT_ME_in_south<0)
		UT_ME_in_south = UT_ME_in_south +360;
	 cosLHA_ME = (sin((-0.583)*pi/180) - sin(35.0*pi/180) * sin(topDecl_ME*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_ME*pi/180));

	if(cosLHA_ME<-1.0){
		cin>>D;
	}
	if(cosLHA_ME>1.0){
		cin>>D;
	}
	 set_ME;
	rise_ME;
	if(cosLHA_ME<1 && cosLHA_ME>-1)
	{
		LHA_ME = acos(cosLHA_ME)*180/pi; //Time in degrees
		if(LHA_ME>(180))
			LHA_ME = LHA_ME- 360;
		if(LHA_ME<-180)
			LHA_ME = LHA_ME + 360;
		set_ME = UT_ME_in_south + LHA_ME;
	}
	while(set_ME>(360))
		set_ME = set_ME - 360;
	while(set_ME<0)
		set_ME = set_ME + 360;
	set_ME = set_ME/15.04107  +4.5;
	perturbed_set_time_ME = set_ME/24;
	}
	perturbed_set_time_ME = perturbed_set_time_ME*24;

	
	//////////////////////////////////////////////////////////////////////////////////PERTURBING RISE

	for(int rise=0;rise<4;rise++){
	d = reald ; //day
	d = d + perturbed_rise_time_M; //perturbing time
	N_ME=48.3313 + 3.24587E-5 * d;  //longitude of the ascending node
    i_ME=7.0047 + 5.00E-8 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    w_ME=29.1241 + 1.01444E-5 * d;  //argument of perihelion
    a_ME=0.387098;  //semi-major axis, or mean distance from Sun
    e_ME=0.205635 + 5.59E-10 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M_ME=168.6562 + 4.0923344368 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_ME>(360))
		M_ME = M_ME - 360;
	while(M_ME<0)
		M_ME = M_ME + 360;
    
	  E_ME = M_ME + e_ME*(180/pi) * sin(M_ME*pi/180) * ( 1.0 + e_ME * cos(M_ME*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_ME>0.05){
        double E0_ME=E_ME;
        double E1_ME=E_ME;
        while(fabs(E0_ME-E1_ME)>0.001){
            E0_ME=E1_ME;
            E1_ME = E_ME - ( E_ME - e_ME*(180/pi) * sin(E_ME*pi/180) - M_ME ) / ( 1 - e_ME * cos(E_ME*pi/180) );
        }
    }
     xv_ME = a_ME * (cos(E_ME*pi/180) - e_ME); //NOT IMPORTANT
     yv_ME = a_ME * sqrt(1.0 - e_ME*e_ME) * sin(E_ME*pi/180); //NOT IMPORTANT
     v_ME = atan2( yv_ME, xv_ME ); //Mercury's true anomaly
     r_ME = sqrt( xv_ME*xv_ME + yv_ME*yv_ME ); //Mercury's distance
     xh_ME = r_ME * ( cos(N_ME*pi/180) * cos(v_ME+w_ME*pi/180) - sin(N_ME*pi/180) * sin(v_ME+w_ME*pi/180) * cos(i_ME*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_ME = r_ME * ( sin(N_ME*pi/180) * cos(v_ME+w_ME*pi/180) + cos(N_ME*pi/180) * sin(v_ME+w_ME*pi/180) * cos(i_ME*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_ME = r_ME * ( sin(v_ME+w_ME*pi/180) * sin(i_ME*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_ME = atan2( yh_ME, xh_ME )*180/pi; //Mercury's ecliptic longitude
     latecl_ME = atan2( zh_ME, sqrt(xh_ME*xh_ME+yh_ME*yh_ME) )*180/pi; //Mercury's ecliptic latitude
	 xg_ME = r_ME * cos(lonecl_ME*pi/180) * cos(latecl_ME*pi/180) +xs; //geocentric position
	 yg_ME = r_ME * sin(lonecl_ME*pi/180) * cos(latecl_ME*pi/180) +ys; //geocentric position
	 zg_ME = r_ME * sin(latecl_ME*pi/180); //geocentric position
	 xe_ME = xg_ME; //equatorial coordinates
     ye_ME = yg_ME * cos(ecl*pi/180) - zg_ME * sin(ecl*pi/180); //equatorial coordinates
     ze_ME = yg_ME * sin(ecl*pi/180) + zg_ME * cos(ecl*pi/180); //equatorial coordinates
     RA_ME  = atan2( ye_ME, xe_ME )*180/pi; //Right Ascension
     Dec_ME = atan2( ze_ME, sqrt(xe_ME*xe_ME+ye_ME*ye_ME) )*180/pi; //Declination
	 par_ME = (8.794/3600) / r_ME;
	 HA_ME = LST - RA_ME; //hour angle of the Moon
	 x_ME = cos(HA_ME*pi/180) * cos(Dec_ME*pi/180);
     y_ME = sin(HA_ME*pi/180) * cos(Dec_ME*pi/180);
     z_ME = sin(Dec_ME*pi/180);
	 xhor_ME = x_ME * sin(35*pi/180) - z_ME * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_ME = y_ME; //Horizontal coordinates
     zhor_ME = x_ME * cos(35*pi/180) + z_ME * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_ME  = atan2( yhor_ME, xhor_ME )*180/pi + 180; //local azimuth
     alt_ME = atan2( zhor_ME, sqrt(xhor_ME*xhor_ME+yhor_ME*yhor_ME) )*180/pi; //local altitude
	 alt_topoc_ME = alt_ME - par_ME * cos(alt_ME*pi/180); //topocentric altitude
	 gclat_ME = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_ME   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_ME = atan( tan(gclat_ME*pi/180) / cos(HA_ME*pi/180) ); //an auxiliary angle
	 topRA_ME   = RA_ME  - par_ME * rho_ME * cos(gclat_ME*pi/180) * sin(HA_ME*pi/180) / cos(Dec_ME*pi/180);
     topDecl_ME = Dec_ME - par_ME * rho_ME * sin(gclat_ME*pi/180) * sin(g_ME - Dec_ME*pi/180) / sin(g_ME);
	
	  UT_ME_in_south = ( topRA_ME - GMST0 - 51.5 );
	while(UT_ME_in_south>360)
		UT_ME_in_south = UT_ME_in_south - 360;
	while(UT_ME_in_south<0)
		UT_ME_in_south = UT_ME_in_south +360;
	 cosLHA_ME = (sin(-0.583*pi/180) - sin(35.0*pi/180) * sin(topDecl_ME*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_ME*pi/180));

	if(cosLHA_ME<-1.0){
		cin>>D;
	}
	if(cosLHA_ME>1.0){
		cin>>D;
	}
	if(cosLHA_ME<1 && cosLHA_ME>-1)
	{
		LHA_ME = acos(cosLHA_ME)*180/pi; //Time in degrees
		if(LHA_ME>(180))
			LHA_ME = LHA_ME- 360;
		if(LHA_ME<-180)
			LHA_ME = LHA_ME + 360;
		rise_ME = UT_ME_in_south - LHA_ME;
	}
	while(rise_ME>(360))
		rise_ME = rise_ME - 360;
	while(rise_ME<0)
		rise_ME = rise_ME + 360;
	rise_ME = rise_ME/15.04107  + 4.5;
	perturbed_rise_time_ME = rise_ME/24;
	}
	perturbed_rise_time_ME = perturbed_rise_time_ME*24;
	UT_ME_in_south = UT_ME_in_south/15;
	UT_ME_in_south = UT_ME_in_south + 4.5;
	double meridianhour_ME = (int)UT_ME_in_south;
	double meridianminute_ME = (int)((UT_ME_in_south - meridianhour_ME) * 60);
	double meridiansecond_ME = (int)((((UT_ME_in_south - meridianhour_ME) * 60) - meridianminute_ME) * 60);
	double sethour_ME = (int)set_ME;
	double setminute_ME = (int)((set_ME - sethour_ME) * 60);
	double setsecond_ME = (int)((((set_ME - sethour_ME) * 60) - setminute_ME) * 60);
	double risehour_ME = (int)rise_ME;
	double riseminute_ME = (int)((rise_ME - risehour_ME) * 60);
	double risesecond_ME = (int)((((rise_ME - risehour_ME) * 60) - riseminute_ME) * 60);
	cout<<"Mercury:"<<endl<< "setting time: "<<sethour_ME<<" : "<<setminute_ME<<" : "<<setsecond_ME<<endl<<"rising time: "<<risehour_ME<<" : "<<riseminute_ME<<" : "<<risesecond_ME<<endl<<"meridian time:"<<meridianhour_ME<<" : "<<meridianminute_ME<<" : "<<meridiansecond_ME<<endl;
		////////////////////////////////Venus's rise and set times////////////////////////////////
    	d = reald + 0.3125;
	N_V=76.6799 + 2.46590E-5 * d;  //longitude of the ascending node
   i_V=3.3946 + 2.75E-8 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
   w_V=54.8910 + 1.38374E-5 * d;  //argument of perihelion
   a_V=0.723330;  //semi-major axis, or mean distance from Sun
   e_V=0.006773 - 1.302E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
   M_V=48.0052 + 1.6021302244 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_V>(360))
		M_V = M_V - 360;
	while(M_V<0)
		M_V = M_V + 360;
    
	  E_V = M_V + e_V*(180/pi) * sin(M_V*pi/180) * ( 1.0 + e_V * cos(M_V*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_V>0.05){
        double E0_V=E_V;
        double E1_V=E_V;
        while(fabs(E0_V-E1_V)>0.001){
            E0_V=E1_V;
            E1_V = E_V - ( E_V - e_V*(180/pi) * sin(E_V*pi/180) - M_V ) / ( 1 - e_V * cos(E_V*pi/180) );
        }
    }
     xv_V = a_V * (cos(E_V*pi/180) - e_V); //NOT IMPORTANT
     yv_V = a_V * sqrt(1.0 - e_V*e_V) * sin(E_V*pi/180); //NOT IMPORTANT
     v_V = atan2( yv_V, xv_V ); //Mercury's true anomaly
     r_V = sqrt( xv_V*xv_V + yv_V*yv_V ); //Mercury's distance
     xh_V = r_V * ( cos(N_V*pi/180) * cos(v_V+w_V*pi/180) - sin(N_V*pi/180) * sin(v_V+w_V*pi/180) * cos(i_V*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_V = r_V * ( sin(N_V*pi/180) * cos(v_V+w_V*pi/180) + cos(N_V*pi/180) * sin(v_V+w_V*pi/180) * cos(i_V*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_V = r_V * ( sin(v_V+w_V*pi/180) * sin(i_V*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_V = atan2( yh_V, xh_V )*180/pi; //Mercury's ecliptic longitude
     latecl_V = atan2( zh_V, sqrt(xh_V*xh_V+yh_V*yh_V) )*180/pi; //Mercury's ecliptic latitude
	 xg_V = r_V * cos(lonecl_V*pi/180) * cos(latecl_V*pi/180) +xs; //geocentric position
	 yg_V = r_V * sin(lonecl_V*pi/180) * cos(latecl_V*pi/180) +ys; //geocentric position
	 zg_V = r_V * sin(latecl_V*pi/180); //geocentric position
	 xe_V = xg_V; //equatorial coordinates
     ye_V = yg_V * cos(ecl*pi/180) - zg_V * sin(ecl*pi/180); //equatorial coordinates
     ze_V = yg_V * sin(ecl*pi/180) + zg_V * cos(ecl*pi/180); //equatorial coordinates
     RA_V  = atan2( ye_V, xe_V )*180/pi; //Right Ascension
     Dec_V = atan2( ze_V, sqrt(xe_V*xe_V+ye_V*ye_V) )*180/pi; //Declination
	 par_V = (8.794/3600) / r_V;
	 HA_V = LST - RA_V; //hour angle of the Moon
	 x_V = cos(HA_V*pi/180) * cos(Dec_V*pi/180);
     y_V = sin(HA_V*pi/180) * cos(Dec_V*pi/180);
     z_V = sin(Dec_V*pi/180);
	 xhor_V = x_V * sin(35*pi/180) - z_V * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_V = y_V; //Horizontal coordinates
     zhor_V = x_V * cos(35*pi/180) + z_V * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_V  = atan2( yhor_V, xhor_V )*180/pi + 180; //local azimuth
     alt_V = atan2( zhor_V, sqrt(xhor_V*xhor_V+yhor_V*yhor_V) )*180/pi; //local altitude
	 alt_topoc_V = alt_V - par_V * cos(alt_V*pi/180); //topocentric altitude
	 gclat_V = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_V   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_V = atan( tan(gclat_V*pi/180) / cos(HA_V*pi/180) ); //an auxiliary angle
	 topRA_V   = RA_V  - par_V * rho_V * cos(gclat_V*pi/180) * sin(HA_V*pi/180) / cos(Dec_V*pi/180);
     topDecl_V = Dec_V - par_V * rho_V * sin(gclat_V*pi/180) * sin(g_V - Dec_V*pi/180) / sin(g_V);
	
	 double UT_V_in_south = ( topRA_V - GMST0 - 51.5 );
	while(UT_V_in_south>360)
		UT_V_in_south = UT_V_in_south - 360;
	while(UT_V_in_south<0)
		UT_V_in_south = UT_V_in_south +360;
	double cosLHA_V = (sin((-0.583)*pi/180) - sin(35.0*pi/180) * sin(topDecl_V*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_V*pi/180));
	double LHA_V;
	if(cosLHA_V<-1.0){
		cin>>D;
	}
	if(cosLHA_V>1.0){
		cin>>D;
	}
	double set_V;
	double rise_V;
	if(cosLHA_V<1 && cosLHA_V>-1)
	{
		LHA_V = acos(cosLHA_V)*180/pi; //Time in degrees
		if(LHA_V>(180))
			LHA_V = LHA_V- 360;
		if(LHA_V<-180)
			LHA_V = LHA_V + 360;
		set_V = UT_V_in_south + LHA_V;
		rise_V = UT_V_in_south - LHA_V;
	}
	while(set_V>(360))
		set_V = set_V - 360;
	while(set_V<0)
		set_V = set_V + 360;
	while(rise_V>(360))
		rise_V = rise_V - 360;
	while(rise_V<0)
		rise_V = rise_V + 360;
	set_V = set_V/15.04107  +4.5;
	rise_V = rise_V/15.04107  +4.5;
	//////////////////////////////////////////////////////////////////////////////////PERTURBING SET
	double perturbed_rise_time_V = rise_V/24;
	double perturbed_set_time_V = set_V/24;

	for(int set=0;set<4;set++){
	d = reald ; //day
	d = d + perturbed_set_time_V; //perturbing time
	N_V=76.6799 + 2.46590E-5 * d;  //longitude of the ascending node
   i_V=3.3946 + 2.75E-8 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
   w_V=54.8910 + 1.38374E-5 * d;  //argument of perihelion
   a_V=0.723330;  //semi-major axis, or mean distance from Sun
   e_V=0.006773 - 1.302E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
   M_V=48.0052 + 1.6021302244 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_V>(360))
		M_V = M_V - 360;
	while(M_V<0)
		M_V = M_V + 360;
    
	  E_V = M_V + e_V*(180/pi) * sin(M_V*pi/180) * ( 1.0 + e_V * cos(M_V*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_V>0.05){
        double E0_V=E_V;
        double E1_V=E_V;
        while(fabs(E0_V-E1_V)>0.001){
            E0_V=E1_V;
            E1_V = E_V - ( E_V - e_V*(180/pi) * sin(E_V*pi/180) - M_V ) / ( 1 - e_V * cos(E_V*pi/180) );
        }
    }
     xv_V = a_V * (cos(E_V*pi/180) - e_V); //NOT IMPORTANT
     yv_V = a_V * sqrt(1.0 - e_V*e_V) * sin(E_V*pi/180); //NOT IMPORTANT
     v_V = atan2( yv_V, xv_V ); //Mercury's true anomaly
     r_V = sqrt( xv_V*xv_V + yv_V*yv_V ); //Mercury's distance
     xh_V = r_V * ( cos(N_V*pi/180) * cos(v_V+w_V*pi/180) - sin(N_V*pi/180) * sin(v_V+w_V*pi/180) * cos(i_V*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_V = r_V * ( sin(N_V*pi/180) * cos(v_V+w_V*pi/180) + cos(N_V*pi/180) * sin(v_V+w_V*pi/180) * cos(i_V*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_V = r_V * ( sin(v_V+w_V*pi/180) * sin(i_V*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_V = atan2( yh_V, xh_V )*180/pi; //Mercury's ecliptic longitude
     latecl_V = atan2( zh_V, sqrt(xh_V*xh_V+yh_V*yh_V) )*180/pi; //Mercury's ecliptic latitude
	 xg_V = r_V * cos(lonecl_V*pi/180) * cos(latecl_V*pi/180) +xs; //geocentric position
	 yg_V = r_V * sin(lonecl_V*pi/180) * cos(latecl_V*pi/180) +ys; //geocentric position
	 zg_V = r_V * sin(latecl_V*pi/180); //geocentric position
	 xe_V = xg_V; //equatorial coordinates
     ye_V = yg_V * cos(ecl*pi/180) - zg_V * sin(ecl*pi/180); //equatorial coordinates
     ze_V = yg_V * sin(ecl*pi/180) + zg_V * cos(ecl*pi/180); //equatorial coordinates
     RA_V  = atan2( ye_V, xe_V )*180/pi; //Right Ascension
     Dec_V = atan2( ze_V, sqrt(xe_V*xe_V+ye_V*ye_V) )*180/pi; //Declination
	 par_V = (8.794/3600) / r_V;
	 HA_V = LST - RA_V; //hour angle of the Moon
	 x_V = cos(HA_V*pi/180) * cos(Dec_V*pi/180);
     y_V = sin(HA_V*pi/180) * cos(Dec_V*pi/180);
     z_V = sin(Dec_V*pi/180);
	 xhor_V = x_V * sin(35*pi/180) - z_V * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_V = y_V; //Horizontal coordinates
     zhor_V = x_V * cos(35*pi/180) + z_V * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_V  = atan2( yhor_V, xhor_V )*180/pi + 180; //local azimuth
     alt_V = atan2( zhor_V, sqrt(xhor_V*xhor_V+yhor_V*yhor_V) )*180/pi; //local altitude
	 alt_topoc_V = alt_V - par_V * cos(alt_V*pi/180); //topocentric altitude
	 gclat_V = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_V   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_V = atan( tan(gclat_V*pi/180) / cos(HA_V*pi/180) ); //an auxiliary angle
	 topRA_V   = RA_V  - par_V * rho_V * cos(gclat_V*pi/180) * sin(HA_V*pi/180) / cos(Dec_V*pi/180);
     topDecl_V = Dec_V - par_V * rho_V * sin(gclat_V*pi/180) * sin(g_V - Dec_V*pi/180) / sin(g_V);
	 
	  UT_V_in_south = ( topRA_V - GMST0 - 51.5 );
	while(UT_V_in_south>360)
		UT_V_in_south = UT_V_in_south - 360;
	while(UT_V_in_south<0)
		UT_V_in_south = UT_V_in_south +360;
	 cosLHA_V = (sin((-0.583)*pi/180) - sin(35.0*pi/180) * sin(topDecl_V*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_V*pi/180));

	if(cosLHA_V<-1.0){
		cin>>D;
	}
	if(cosLHA_V>1.0){
		cin>>D;
	}
	 set_V;
	rise_V;
	if(cosLHA_V<1 && cosLHA_V>-1)
	{
		LHA_V = acos(cosLHA_V)*180/pi; //Time in degrees
		if(LHA_V>(180))
			LHA_V = LHA_V- 360;
		if(LHA_V<-180)
			LHA_V = LHA_V + 360;
		set_V = UT_V_in_south + LHA_V;
	}
	while(set_V>(360))
		set_V = set_V - 360;
	while(set_V<0)
		set_V = set_V + 360;
	set_V = set_V/15.04107  +4.5;
	perturbed_set_time_V = set_V/24;
	}
	perturbed_set_time_V = perturbed_set_time_V*24;

	
	//////////////////////////////////////////////////////////////////////////////////PERTURBING RISE

	for(int rise=0;rise<4;rise++){
	d = reald ; //day
	d = d + perturbed_rise_time_M; //perturbing time
	N_V=76.6799 + 2.46590E-5 * d;  //longitude of the ascending node
   i_V=3.3946 + 2.75E-8 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
   w_V=54.8910 + 1.38374E-5 * d;  //argument of perihelion
   a_V=0.723330;  //semi-major axis, or mean distance from Sun
   e_V=0.006773 - 1.302E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
   M_V=48.0052 + 1.6021302244 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_V>(360))
		M_V = M_V - 360;
	while(M_V<0)
		M_V = M_V + 360;
    
	  E_V = M_V + e_V*(180/pi) * sin(M_V*pi/180) * ( 1.0 + e_V * cos(M_V*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_V>0.05){
        double E0_V=E_V;
        double E1_V=E_V;
        while(fabs(E0_V-E1_V)>0.001){
            E0_V=E1_V;
            E1_V = E_V - ( E_V - e_V*(180/pi) * sin(E_V*pi/180) - M_V ) / ( 1 - e_V * cos(E_V*pi/180) );
        }
    }
     xv_V = a_V * (cos(E_V*pi/180) - e_V); //NOT IMPORTANT
     yv_V = a_V * sqrt(1.0 - e_V*e_V) * sin(E_V*pi/180); //NOT IMPORTANT
     v_V = atan2( yv_V, xv_V ); //Mercury's true anomaly
     r_V = sqrt( xv_V*xv_V + yv_V*yv_V ); //Mercury's distance
     xh_V = r_V * ( cos(N_V*pi/180) * cos(v_V+w_V*pi/180) - sin(N_V*pi/180) * sin(v_V+w_V*pi/180) * cos(i_V*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_V = r_V * ( sin(N_V*pi/180) * cos(v_V+w_V*pi/180) + cos(N_V*pi/180) * sin(v_V+w_V*pi/180) * cos(i_V*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_V = r_V * ( sin(v_V+w_V*pi/180) * sin(i_V*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_V = atan2( yh_V, xh_V )*180/pi; //Mercury's ecliptic longitude
     latecl_V = atan2( zh_V, sqrt(xh_V*xh_V+yh_V*yh_V) )*180/pi; //Mercury's ecliptic latitude
	 xg_V = r_V * cos(lonecl_V*pi/180) * cos(latecl_V*pi/180) +xs; //geocentric position
	 yg_V = r_V * sin(lonecl_V*pi/180) * cos(latecl_V*pi/180) +ys; //geocentric position
	 zg_V = r_V * sin(latecl_V*pi/180); //geocentric position
	 xe_V = xg_V; //equatorial coordinates
     ye_V = yg_V * cos(ecl*pi/180) - zg_V * sin(ecl*pi/180); //equatorial coordinates
     ze_V = yg_V * sin(ecl*pi/180) + zg_V * cos(ecl*pi/180); //equatorial coordinates
     RA_V  = atan2( ye_V, xe_V )*180/pi; //Right Ascension
     Dec_V = atan2( ze_V, sqrt(xe_V*xe_V+ye_V*ye_V) )*180/pi; //Declination
	 par_V = (8.794/3600) / r_V;
	 HA_V = LST - RA_V; //hour angle of the Moon
	 x_V = cos(HA_V*pi/180) * cos(Dec_V*pi/180);
     y_V = sin(HA_V*pi/180) * cos(Dec_V*pi/180);
     z_V = sin(Dec_V*pi/180);
	 xhor_V = x_V * sin(35*pi/180) - z_V * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_V = y_V; //Horizontal coordinates
     zhor_V = x_V * cos(35*pi/180) + z_V * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_V  = atan2( yhor_V, xhor_V )*180/pi + 180; //local azimuth
     alt_V = atan2( zhor_V, sqrt(xhor_V*xhor_V+yhor_V*yhor_V) )*180/pi; //local altitude
	 alt_topoc_V = alt_V - par_V * cos(alt_V*pi/180); //topocentric altitude
	 gclat_V = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_V   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_V = atan( tan(gclat_V*pi/180) / cos(HA_V*pi/180) ); //an auxiliary angle
	 topRA_V   = RA_V  - par_V * rho_V * cos(gclat_V*pi/180) * sin(HA_V*pi/180) / cos(Dec_V*pi/180);
     topDecl_V = Dec_V - par_V * rho_V * sin(gclat_V*pi/180) * sin(g_V - Dec_V*pi/180) / sin(g_V);
	
	  UT_V_in_south = ( topRA_V - GMST0 - 51.5 );
	while(UT_V_in_south>360)
		UT_V_in_south = UT_V_in_south - 360;
	while(UT_V_in_south<0)
		UT_V_in_south = UT_V_in_south +360;
	 cosLHA_V = (sin(-0.583*pi/180) - sin(35.0*pi/180) * sin(topDecl_V*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_V*pi/180));

	if(cosLHA_V<-1.0){
		cin>>D;
	}
	if(cosLHA_V>1.0){
		cin>>D;
	}
	if(cosLHA_V<1 && cosLHA_V>-1)
	{
		LHA_V = acos(cosLHA_V)*180/pi; //Time in degrees
		if(LHA_V>(180))
			LHA_V = LHA_V- 360;
		if(LHA_V<-180)
			LHA_V = LHA_V + 360;
		rise_V = UT_V_in_south - LHA_V;
	}
	while(rise_V>(360))
		rise_V = rise_V - 360;
	while(rise_V<0)
		rise_V = rise_V + 360;
	rise_V = rise_V/15.04107  + 4.5;
	perturbed_rise_time_V = rise_V/24;
	}
	perturbed_rise_time_V = perturbed_rise_time_V*24;
	UT_V_in_south = UT_V_in_south/15;
	UT_V_in_south = UT_V_in_south + 4.5;
	double meridianhour_V = (int)UT_V_in_south;
	double meridianminute_V = (int)((UT_V_in_south - meridianhour_V) * 60);
	double meridiansecond_V = (int)((((UT_V_in_south - meridianhour_V) * 60) - meridianminute_V) * 60);
	double sethour_V = (int)set_V;
	double setminute_V = (int)((set_V - sethour_V) * 60);
	double setsecond_V = (int)((((set_V - sethour_V) * 60) - setminute_V) * 60);
	double risehour_V = (int)rise_V;
	double riseminute_V = (int)((rise_V - risehour_V) * 60);
	double risesecond_V = (int)((((rise_V - risehour_V) * 60) - riseminute_V) * 60);
	cout<<"Venus:"<<endl<< "setting time: "<<sethour_V<<" : "<<setminute_V<<" : "<<setsecond_V<<endl<<"rising time: "<<risehour_V<<" : "<<riseminute_V<<" : "<<risesecond_V<<endl<<"meridian time:"<<meridianhour_V<<" : "<<meridianminute_V<<" : "<<meridiansecond_V<<endl;

	////////////////////////////////Mars's rise and set times////////////////////////////////
    d = reald + 0.3125;
	N_MA=49.5574 + 2.11081E-5 * d;  //longitude of the ascending node
    i_MA=1.8497 - 1.78E-8 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    w_MA=286.5016 + 2.92961E-5 * d;  //argument of perihelion
    a_MA=1.523688;  //semi-major axis, or mean distance from Sun
    e_MA=0.093405 + 2.516E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M_MA=18.6021 + 0.5240207766 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_MA>(360))
		M_MA = M_MA - 360;
	while(M_MA<0)
		M_MA = M_MA + 360;
   	
	  E_MA = M_MA + e_MA*(180/pi) * sin(M_MA*pi/180) * ( 1.0 + e_MA * cos(M_MA*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_MA>0.05){
        double E0_MA=E_MA;
        double E1_MA=E_MA;
        while(fabs(E0_MA-E1_MA)>0.001){
            E0_MA=E1_MA;
            E1_MA = E_MA - ( E_MA - e_MA*(180/pi) * sin(E_MA*pi/180) - M_MA ) / ( 1 - e_MA * cos(E_MA*pi/180) );
        }
    }
     xv_MA = a_MA * (cos(E_MA*pi/180) - e_MA); //NOT IMPORTANT
     yv_MA = a_MA * sqrt(1.0 - e_MA*e_MA) * sin(E_MA*pi/180); //NOT IMPORTANT
     v_MA = atan2( yv_MA, xv_MA ); //Mercury's true anomaly
     r_MA = sqrt( xv_MA*xv_MA + yv_MA*yv_MA ); //Mercury's distance
     xh_MA = r_MA * ( cos(N_MA*pi/180) * cos(v_MA+w_MA*pi/180) - sin(N_MA*pi/180) * sin(v_MA+w_MA*pi/180) * cos(i_MA*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_MA = r_MA * ( sin(N_MA*pi/180) * cos(v_MA+w_MA*pi/180) + cos(N_MA*pi/180) * sin(v_MA+w_MA*pi/180) * cos(i_MA*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_MA = r_MA * ( sin(v_MA+w_MA*pi/180) * sin(i_MA*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_MA = atan2( yh_MA, xh_MA )*180/pi; //Mercury's ecliptic longitude
     latecl_MA = atan2( zh_MA, sqrt(xh_MA*xh_MA+yh_MA*yh_MA) )*180/pi; //Mercury's ecliptic latitude
	 xg_MA = r_MA * cos(lonecl_MA*pi/180) * cos(latecl_MA*pi/180) +xs; //geocentric position
	 yg_MA = r_MA * sin(lonecl_MA*pi/180) * cos(latecl_MA*pi/180) +ys; //geocentric position
	 zg_MA = r_MA * sin(latecl_MA*pi/180); //geocentric position
	 xe_MA = xg_MA; //equatorial coordinates
     ye_MA = yg_MA * cos(ecl*pi/180) - zg_MA * sin(ecl*pi/180); //equatorial coordinates
     ze_MA = yg_MA * sin(ecl*pi/180) + zg_MA * cos(ecl*pi/180); //equatorial coordinates
     RA_MA  = atan2( ye_MA, xe_MA )*180/pi; //Right Ascension
     Dec_MA = atan2( ze_MA, sqrt(xe_MA*xe_MA+ye_MA*ye_MA) )*180/pi; //Declination
	 par_MA = (8.794/3600) / r_MA;
	 HA_MA = LST - RA_MA; //hour angle of the Moon
	 x_MA = cos(HA_MA*pi/180) * cos(Dec_MA*pi/180);
     y_MA = sin(HA_MA*pi/180) * cos(Dec_MA*pi/180);
     z_MA = sin(Dec_MA*pi/180);
	 xhor_MA = x_MA * sin(35*pi/180) - z_MA * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_MA = y_MA; //Horizontal coordinates
     zhor_MA = x_MA * cos(35*pi/180) + z_MA * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_MA  = atan2( yhor_MA, xhor_MA )*180/pi + 180; //local azimuth
     alt_MA = atan2( zhor_MA, sqrt(xhor_MA*xhor_MA+yhor_MA*yhor_MA) )*180/pi; //local altitude
	 alt_topoc_MA = alt_MA - par_MA * cos(alt_MA*pi/180); //topocentric altitude
	 gclat_MA = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_MA   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_MA = atan( tan(gclat_MA*pi/180) / cos(HA_MA*pi/180) ); //an auxiliary angle
	 topRA_MA   = RA_MA  - par_MA * rho_MA * cos(gclat_MA*pi/180) * sin(HA_MA*pi/180) / cos(Dec_MA*pi/180);
     topDecl_MA = Dec_MA - par_MA * rho_MA * sin(gclat_MA*pi/180) * sin(g_MA - Dec_MA*pi/180) / sin(g_MA);
	
	 double UT_MA_in_south = ( topRA_MA - GMST0 - 51.5 );
	while(UT_MA_in_south>360)
		UT_MA_in_south = UT_MA_in_south - 360;
	while(UT_MA_in_south<0)
		UT_MA_in_south = UT_MA_in_south +360;
	double cosLHA_MA = (sin((-0.583)*pi/180) - sin(35.0*pi/180) * sin(topDecl_MA*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_MA*pi/180));
	double LHA_MA;
	if(cosLHA_MA<-1.0){
		cin>>D;
	}
	if(cosLHA_MA>1.0){
		cin>>D;
	}
	double set_MA;
	double rise_MA;
	if(cosLHA_MA<1 && cosLHA_MA>-1)
	{
		LHA_MA = acos(cosLHA_MA)*180/pi; //Time in degrees
		if(LHA_MA>(180))
			LHA_MA = LHA_MA- 360;
		if(LHA_MA<-180)
			LHA_MA = LHA_MA + 360;
		set_MA = UT_MA_in_south + LHA_MA;
		rise_MA = UT_MA_in_south - LHA_MA;
	}
	while(set_MA>(360))
		set_MA = set_MA - 360;
	while(set_MA<0)
		set_MA = set_MA + 360;
	while(rise_MA>(360))
		rise_MA = rise_MA - 360;
	while(rise_MA<0)
		rise_MA = rise_MA + 360;
	set_MA = set_MA/15.04107  +4.5;
	rise_MA = rise_MA/15.04107  +4.5;
	//////////////////////////////////////////////////////////////////////////////////PERTURBING SET
	double perturbed_rise_time_MA = rise_MA/24;
	double perturbed_set_time_MA = set_MA/24;

	for(int set=0;set<4;set++){
	d = reald ; //day
	d = d + perturbed_set_time_MA; //perturbing time
	N_MA=49.5574 + 2.11081E-5 * d;  //longitude of the ascending node
    i_MA=1.8497 - 1.78E-8 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    w_MA=286.5016 + 2.92961E-5 * d;  //argument of perihelion
    a_MA=1.523688;  //semi-major axis, or mean distance from Sun
    e_MA=0.093405 + 2.516E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M_MA=18.6021 + 0.5240207766 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_MA>(360))
		M_MA = M_MA - 360;
	while(M_MA<0)
		M_MA = M_MA + 360;
   	
	  E_MA = M_MA + e_MA*(180/pi) * sin(M_MA*pi/180) * ( 1.0 + e_MA * cos(M_MA*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_MA>0.05){
        double E0_MA=E_MA;
        double E1_MA=E_MA;
        while(fabs(E0_MA-E1_MA)>0.001){
            E0_MA=E1_MA;
            E1_MA = E_MA - ( E_MA - e_MA*(180/pi) * sin(E_MA*pi/180) - M_MA ) / ( 1 - e_MA * cos(E_MA*pi/180) );
        }
    }
     xv_MA = a_MA * (cos(E_MA*pi/180) - e_MA); //NOT IMPORTANT
     yv_MA = a_MA * sqrt(1.0 - e_MA*e_MA) * sin(E_MA*pi/180); //NOT IMPORTANT
     v_MA = atan2( yv_MA, xv_MA ); //Mercury's true anomaly
     r_MA = sqrt( xv_MA*xv_MA + yv_MA*yv_MA ); //Mercury's distance
     xh_MA = r_MA * ( cos(N_MA*pi/180) * cos(v_MA+w_MA*pi/180) - sin(N_MA*pi/180) * sin(v_MA+w_MA*pi/180) * cos(i_MA*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_MA = r_MA * ( sin(N_MA*pi/180) * cos(v_MA+w_MA*pi/180) + cos(N_MA*pi/180) * sin(v_MA+w_MA*pi/180) * cos(i_MA*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_MA = r_MA * ( sin(v_MA+w_MA*pi/180) * sin(i_MA*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_MA = atan2( yh_MA, xh_MA )*180/pi; //Mercury's ecliptic longitude
     latecl_MA = atan2( zh_MA, sqrt(xh_MA*xh_MA+yh_MA*yh_MA) )*180/pi; //Mercury's ecliptic latitude
	 xg_MA = r_MA * cos(lonecl_MA*pi/180) * cos(latecl_MA*pi/180) +xs; //geocentric position
	 yg_MA = r_MA * sin(lonecl_MA*pi/180) * cos(latecl_MA*pi/180) +ys; //geocentric position
	 zg_MA = r_MA * sin(latecl_MA*pi/180); //geocentric position
	 xe_MA = xg_MA; //equatorial coordinates
     ye_MA = yg_MA * cos(ecl*pi/180) - zg_MA * sin(ecl*pi/180); //equatorial coordinates
     ze_MA = yg_MA * sin(ecl*pi/180) + zg_MA * cos(ecl*pi/180); //equatorial coordinates
     RA_MA  = atan2( ye_MA, xe_MA )*180/pi; //Right Ascension
     Dec_MA = atan2( ze_MA, sqrt(xe_MA*xe_MA+ye_MA*ye_MA) )*180/pi; //Declination
	 par_MA = (8.794/3600) / r_MA;
	 HA_MA = LST - RA_MA; //hour angle of the Moon
	 x_MA = cos(HA_MA*pi/180) * cos(Dec_MA*pi/180);
     y_MA = sin(HA_MA*pi/180) * cos(Dec_MA*pi/180);
     z_MA = sin(Dec_MA*pi/180);
	 xhor_MA = x_MA * sin(35*pi/180) - z_MA * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_MA = y_MA; //Horizontal coordinates
     zhor_MA = x_MA * cos(35*pi/180) + z_MA * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_MA  = atan2( yhor_MA, xhor_MA )*180/pi + 180; //local azimuth
     alt_MA = atan2( zhor_MA, sqrt(xhor_MA*xhor_MA+yhor_MA*yhor_MA) )*180/pi; //local altitude
	 alt_topoc_MA = alt_MA - par_MA * cos(alt_MA*pi/180); //topocentric altitude
	 gclat_MA = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_MA   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_MA = atan( tan(gclat_MA*pi/180) / cos(HA_MA*pi/180) ); //an auxiliary angle
	 topRA_MA   = RA_MA  - par_MA * rho_MA * cos(gclat_MA*pi/180) * sin(HA_MA*pi/180) / cos(Dec_MA*pi/180);
     topDecl_MA = Dec_MA - par_MA * rho_MA * sin(gclat_MA*pi/180) * sin(g_MA - Dec_MA*pi/180) / sin(g_MA);
	 
	  UT_MA_in_south = ( topRA_MA - GMST0 - 51.5 );
	while(UT_MA_in_south>360)
		UT_MA_in_south = UT_MA_in_south - 360;
	while(UT_MA_in_south<0)
		UT_MA_in_south = UT_MA_in_south +360;
	 cosLHA_MA = (sin((-0.583)*pi/180) - sin(35.0*pi/180) * sin(topDecl_MA*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_MA*pi/180));

	if(cosLHA_MA<-1.0){
		cin>>D;
	}
	if(cosLHA_MA>1.0){
		cin>>D;
	}
	 set_MA;
	rise_MA;
	if(cosLHA_MA<1 && cosLHA_MA>-1)
	{
		LHA_MA = acos(cosLHA_MA)*180/pi; //Time in degrees
		if(LHA_MA>(180))
			LHA_MA = LHA_MA- 360;
		if(LHA_MA<-180)
			LHA_MA = LHA_MA + 360;
		set_MA = UT_MA_in_south + LHA_MA;
	}
	while(set_MA>(360))
		set_MA = set_MA - 360;
	while(set_MA<0)
		set_MA = set_MA + 360;
	set_MA = set_MA/15.04107  +4.5;
	perturbed_set_time_MA = set_MA/24;
	}
	perturbed_set_time_MA = perturbed_set_time_MA*24;

	
	//////////////////////////////////////////////////////////////////////////////////PERTURBING RISE

	for(int rise=0;rise<4;rise++){
	d = reald ; //day
	d = d + perturbed_rise_time_M; //perturbing time
	N_MA=49.5574 + 2.11081E-5 * d;  //longitude of the ascending node
    i_MA=1.8497 - 1.78E-8 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    w_MA=286.5016 + 2.92961E-5 * d;  //argument of perihelion
    a_MA=1.523688;  //semi-major axis, or mean distance from Sun
    e_MA=0.093405 + 2.516E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M_MA=18.6021 + 0.5240207766 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_MA>(360))
		M_MA = M_MA - 360;
	while(M_MA<0)
		M_MA = M_MA + 360;
   		while(M_MA>(360))
		M_MA = M_MA - 360;
	while(M_MA<0)
		M_MA = M_MA + 360;
    
	  E_MA = M_MA + e_MA*(180/pi) * sin(M_MA*pi/180) * ( 1.0 + e_MA * cos(M_MA*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_MA>0.05){
        double E0_MA=E_MA;
        double E1_MA=E_MA;
        while(fabs(E0_MA-E1_MA)>0.001){
            E0_MA=E1_MA;
            E1_MA = E_MA - ( E_MA - e_MA*(180/pi) * sin(E_MA*pi/180) - M_MA ) / ( 1 - e_MA * cos(E_MA*pi/180) );
        }
    }
     xv_MA = a_MA * (cos(E_MA*pi/180) - e_MA); //NOT IMPORTANT
     yv_MA = a_MA * sqrt(1.0 - e_MA*e_MA) * sin(E_MA*pi/180); //NOT IMPORTANT
     v_MA = atan2( yv_MA, xv_MA ); //Mercury's true anomaly
     r_MA = sqrt( xv_MA*xv_MA + yv_MA*yv_MA ); //Mercury's distance
     xh_MA = r_MA * ( cos(N_MA*pi/180) * cos(v_MA+w_MA*pi/180) - sin(N_MA*pi/180) * sin(v_MA+w_MA*pi/180) * cos(i_MA*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_MA = r_MA * ( sin(N_MA*pi/180) * cos(v_MA+w_MA*pi/180) + cos(N_MA*pi/180) * sin(v_MA+w_MA*pi/180) * cos(i_MA*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_MA = r_MA * ( sin(v_MA+w_MA*pi/180) * sin(i_MA*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_MA = atan2( yh_MA, xh_MA )*180/pi; //Mercury's ecliptic longitude
     latecl_MA = atan2( zh_MA, sqrt(xh_MA*xh_MA+yh_MA*yh_MA) )*180/pi; //Mercury's ecliptic latitude
	 xg_MA = r_MA * cos(lonecl_MA*pi/180) * cos(latecl_MA*pi/180) +xs; //geocentric position
	 yg_MA = r_MA * sin(lonecl_MA*pi/180) * cos(latecl_MA*pi/180) +ys; //geocentric position
	 zg_MA = r_MA * sin(latecl_MA*pi/180); //geocentric position
	 xe_MA = xg_MA; //equatorial coordinates
     ye_MA = yg_MA * cos(ecl*pi/180) - zg_MA * sin(ecl*pi/180); //equatorial coordinates
     ze_MA = yg_MA * sin(ecl*pi/180) + zg_MA * cos(ecl*pi/180); //equatorial coordinates
     RA_MA  = atan2( ye_MA, xe_MA )*180/pi; //Right Ascension
     Dec_MA = atan2( ze_MA, sqrt(xe_MA*xe_MA+ye_MA*ye_MA) )*180/pi; //Declination
	 par_MA = (8.794/3600) / r_MA;
	 HA_MA = LST - RA_MA; //hour angle of the Moon
	 x_MA = cos(HA_MA*pi/180) * cos(Dec_MA*pi/180);
     y_MA = sin(HA_MA*pi/180) * cos(Dec_MA*pi/180);
     z_MA = sin(Dec_MA*pi/180);
	 xhor_MA = x_MA * sin(35*pi/180) - z_MA * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_MA = y_MA; //Horizontal coordinates
     zhor_MA = x_MA * cos(35*pi/180) + z_MA * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_MA  = atan2( yhor_MA, xhor_MA )*180/pi + 180; //local azimuth
     alt_MA = atan2( zhor_MA, sqrt(xhor_MA*xhor_MA+yhor_MA*yhor_MA) )*180/pi; //local altitude
	 alt_topoc_MA = alt_MA - par_MA * cos(alt_MA*pi/180); //topocentric altitude
	 gclat_MA = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_MA   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_MA = atan( tan(gclat_MA*pi/180) / cos(HA_MA*pi/180) ); //an auxiliary angle
	 topRA_MA   = RA_MA  - par_MA * rho_MA * cos(gclat_MA*pi/180) * sin(HA_MA*pi/180) / cos(Dec_MA*pi/180);
     topDecl_MA = Dec_MA - par_MA * rho_MA * sin(gclat_MA*pi/180) * sin(g_MA - Dec_MA*pi/180) / sin(g_MA);
	
	  UT_MA_in_south = ( topRA_MA - GMST0 - 51.5 );
	while(UT_MA_in_south>360)
		UT_MA_in_south = UT_MA_in_south - 360;
	while(UT_MA_in_south<0)
		UT_MA_in_south = UT_MA_in_south +360;
	 cosLHA_MA = (sin(-0.583*pi/180) - sin(35.0*pi/180) * sin(topDecl_MA*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_MA*pi/180));

	if(cosLHA_MA<-1.0){
		cin>>D;
	}
	if(cosLHA_MA>1.0){
		cin>>D;
	}
	if(cosLHA_MA<1 && cosLHA_MA>-1)
	{
		LHA_MA = acos(cosLHA_MA)*180/pi; //Time in degrees
		if(LHA_MA>(180))
			LHA_MA = LHA_MA- 360;
		if(LHA_MA<-180)
			LHA_MA = LHA_MA + 360;
		rise_MA = UT_MA_in_south - LHA_MA;
	}
	while(rise_MA>(360))
		rise_MA = rise_MA - 360;
	while(rise_MA<0)
		rise_MA = rise_MA + 360;
	rise_MA = rise_MA/15.04107  + 4.5;
	perturbed_rise_time_MA = rise_MA/24;
	}
	perturbed_rise_time_MA = perturbed_rise_time_MA*24;
	UT_MA_in_south = UT_MA_in_south/15;
	UT_MA_in_south = UT_MA_in_south + 4.5;
	double meridianhour_MA = (int)UT_MA_in_south;
	double meridianminute_MA = (int)((UT_MA_in_south - meridianhour_MA) * 60);
	double meridiansecond_MA = (int)((((UT_MA_in_south - meridianhour_MA) * 60) - meridianminute_MA) * 60);
	double sethour_MA = (int)set_MA;
	double setminute_MA = (int)((set_MA - sethour_MA) * 60);
	double setsecond_MA = (int)((((set_MA - sethour_MA) * 60) - setminute_MA) * 60);
	double risehour_MA = (int)rise_MA;
	double riseminute_MA = (int)((rise_MA - risehour_MA) * 60);
	double risesecond_MA = (int)((((rise_MA - risehour_MA) * 60) - riseminute_MA) * 60);
	cout<<"Mars:"<<endl<< "setting time: "<<sethour_MA<<" : "<<setminute_MA<<" : "<<setsecond_MA<<endl<<"rising time: "<<risehour_MA<<" : "<<riseminute_MA<<" : "<<risesecond_MA<<endl<<"meridian time:"<<meridianhour_MA<<" : "<<meridianminute_MA<<" : "<<meridiansecond_MA<<endl;

	////////////////////////////////Jupiter's rise and set times////////////////////////////////
    d = reald + 0.3125;
	N_J=100.4542 + 2.76854E-5 * d;  //longitude of the ascending node
    i_J=1.3030 - 1.557E-7 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    w_J=273.8777 + 1.64505E-5 * d;  //argument of perihelion
    a_J=5.20256  ;  //semi-major axis, or mean distance from Sun
    e_J=0.048498 + 4.469E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M_J=19.8950 + 0.0830853001 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_J>(360))
		M_J = M_J - 360;
	while(M_J<0)
		M_J = M_J + 360;
    
	  E_J = M_J + e_J*(180/pi) * sin(M_J*pi/180) * ( 1.0 + e_J * cos(M_J*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_J>0.05){
        double E0_J=E_J;
        double E1_J=E_J;
        while(fabs(E0_J-E1_J)>0.001){
            E0_J=E1_J;
            E1_J = E_J - ( E_J - e_J*(180/pi) * sin(E_J*pi/180) - M_J ) / ( 1 - e_J * cos(E_J*pi/180) );
        }
    }
     xv_J = a_J * (cos(E_J*pi/180) - e_J); //NOT IMPORTANT
     yv_J = a_J * sqrt(1.0 - e_J*e_J) * sin(E_J*pi/180); //NOT IMPORTANT
     v_J = atan2( yv_J, xv_J ); //Mercury's true anomaly
     r_J = sqrt( xv_J*xv_J + yv_J*yv_J ); //Mercury's distance
     xh_J = r_J * ( cos(N_J*pi/180) * cos(v_J+w_J*pi/180) - sin(N_J*pi/180) * sin(v_J+w_J*pi/180) * cos(i_J*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_J = r_J * ( sin(N_J*pi/180) * cos(v_J+w_J*pi/180) + cos(N_J*pi/180) * sin(v_J+w_J*pi/180) * cos(i_J*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_J = r_J * ( sin(v_J+w_J*pi/180) * sin(i_J*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_J = atan2( yh_J, xh_J )*180/pi; //Mercury's ecliptic longitude
     latecl_J = atan2( zh_J, sqrt(xh_J*xh_J+yh_J*yh_J) )*180/pi; //Mercury's ecliptic latitude
	 M_SA=316.9670 + 0.0334442282 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_S>(360))
		M_S = M_S - 360;
	while(M_S<0)
		M_S = M_S + 360;
	 M_U=142.5905 + 0.011725806 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_U>(360))
		M_U = M_U - 360;
	while(M_U<0)
		M_U = M_U + 360;
	
    
	 lonecl_J = lonecl_J -0.332 * sin(2*M_J*pi/180 - 5*M_SA*pi/180 - 67.6*pi/180 )
    -0.056 * sin(2*M_J*pi/180 - 2*M_SA*pi/180 + 21*pi/180 )
    +0.042 * sin(3*M_J*pi/180 - 5*M_SA*pi/180 + 21*pi/180 )
    -0.036 * sin(M_J*pi/180 - 2*M_SA*pi/180)
    +0.022 * cos(M_J*pi/180 - M_SA*pi/180)
    +0.023 * sin(2*M_J*pi/180 - 3*M_SA + 52*pi/180 )
    -0.016 * sin(M_J*pi/180 - 5*M_SA*pi/180 - 69*pi/180 );

	lonecl_SA = lonecl_SA +0.812 * sin(2*M_J*pi/180 - 5*M_SA*pi/180 - 67.6*pi/180 )
    -0.229 * cos(2*M_J*pi/180 - 4*M_SA*pi/180 - 2*pi/180 )
    +0.119 * sin(M_J*pi/180 - 2*M_SA*pi/180 - 3*pi/180 )
    +0.046 * sin(2*M_J*pi/180 - 6*M_SA*pi/180 - 69*pi/180 )
    +0.014 * sin(M_J*pi/180 - 3*M_SA*pi/180 + 32*pi/180 );

	latecl_SA = latecl_SA -0.020 * cos(2*M_J*pi/180 - 4*M_SA*pi/180 - 2*pi/180 )
    +0.018 * sin(2*M_J*pi/180 - 6*M_SA*pi/180 - 49*pi/180 );

	lonecl_U = lonecl_U +0.040 * sin(M_SA*pi/180 - 2*M_U*pi/180 + 6*pi/180 )
    +0.035 * sin(M_SA*pi/180 - 3*M_U*pi/180 + 33*pi/180 )
    -0.015 * sin(M_J*pi/180 - M_U*pi/180 + 20*pi/180 );

	 xg_J = r_J * cos(lonecl_J*pi/180) * cos(latecl_J*pi/180) +xs; //geocentric position
	 yg_J = r_J * sin(lonecl_J*pi/180) * cos(latecl_J*pi/180) +ys; //geocentric position
	 zg_J = r_J * sin(latecl_J*pi/180); //geocentric position
	 xe_J = xg_J; //equatorial coordinates
     ye_J = yg_J * cos(ecl*pi/180) - zg_J * sin(ecl*pi/180); //equatorial coordinates
     ze_J = yg_J * sin(ecl*pi/180) + zg_J * cos(ecl*pi/180); //equatorial coordinates
     RA_J  = atan2( ye_J, xe_J )*180/pi; //Right Ascension
     Dec_J = atan2( ze_J, sqrt(xe_J*xe_J+ye_J*ye_J) )*180/pi; //Declination
	 par_J = (8.794/3600) / r_J;
	 HA_J = LST - RA_J; //hour angle of the Moon
	 x_J = cos(HA_J*pi/180) * cos(Dec_J*pi/180);
     y_J = sin(HA_J*pi/180) * cos(Dec_J*pi/180);
     z_J = sin(Dec_J*pi/180);
	 xhor_J = x_J * sin(35*pi/180) - z_J * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_J = y_J; //Horizontal coordinates
     zhor_J = x_J * cos(35*pi/180) + z_J * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_J  = atan2( yhor_J, xhor_J )*180/pi + 180; //local azimuth
     alt_J = atan2( zhor_J, sqrt(xhor_J*xhor_J+yhor_J*yhor_J) )*180/pi; //local altitude
	 alt_topoc_J = alt_J - par_J * cos(alt_J*pi/180); //topocentric altitude
	 gclat_J = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_J   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_J = atan( tan(gclat_J*pi/180) / cos(HA_J*pi/180) ); //an auxiliary angle
	 topRA_J   = RA_J  - par_J * rho_J * cos(gclat_J*pi/180) * sin(HA_J*pi/180) / cos(Dec_J*pi/180);
     topDecl_J = Dec_J - par_J * rho_J * sin(gclat_J*pi/180) * sin(g_J - Dec_J*pi/180) / sin(g_J);
	
	 double UT_J_in_south = ( topRA_J - GMST0 - 51.5 );
	while(UT_J_in_south>360)
		UT_J_in_south = UT_J_in_south - 360;
	while(UT_J_in_south<0)
		UT_J_in_south = UT_J_in_south +360;
	double cosLHA_J = (sin((-0.583)*pi/180) - sin(35.0*pi/180) * sin(topDecl_J*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_J*pi/180));
	double LHA_J;
	if(cosLHA_J<-1.0){
		cin>>D;
	}
	if(cosLHA_J>1.0){
		cin>>D;
	}
	double set_J;
	double rise_J;
	if(cosLHA_J<1 && cosLHA_J>-1)
	{
		LHA_J = acos(cosLHA_J)*180/pi; //Time in degrees
		if(LHA_J>(180))
			LHA_J = LHA_J- 360;
		if(LHA_J<-180)
			LHA_J = LHA_J + 360;
		set_J = UT_J_in_south + LHA_J;
		rise_J = UT_J_in_south - LHA_J;
	}
	while(set_J>(360))
		set_J = set_J - 360;
	while(set_J<0)
		set_J = set_J + 360;
	while(rise_J>(360))
		rise_J = rise_J - 360;
	while(rise_J<0)
		rise_J = rise_J + 360;
	set_J = set_J/15.04107  +4.5;
	rise_J = rise_J/15.04107  +4.5;
	//////////////////////////////////////////////////////////////////////////////////PERTURBING SET
	double perturbed_rise_time_J = rise_J/24;
	double perturbed_set_time_J = set_J/24;

	for(int set=0;set<4;set++){
	d = reald ; //day
	d = d + perturbed_set_time_J; //perturbing time
	N_J=100.4542 + 2.76854E-5 * d;  //longitude of the ascending node
    i_J=1.3030 - 1.557E-7 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    w_J=273.8777 + 1.64505E-5 * d;  //argument of perihelion
    a_J=5.20256  ;  //semi-major axis, or mean distance from Sun
    e_J=0.048498 + 4.469E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M_J=19.8950 + 0.0830853001 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_J>(360))
		M_J = M_J - 360;
	while(M_J<0)
		M_J = M_J + 360;
    
	  E_J = M_J + e_J*(180/pi) * sin(M_J*pi/180) * ( 1.0 + e_J * cos(M_J*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_J>0.05){
        double E0_J=E_J;
        double E1_J=E_J;
        while(fabs(E0_J-E1_J)>0.001){
            E0_J=E1_J;
            E1_J = E_J - ( E_J - e_J*(180/pi) * sin(E_J*pi/180) - M_J ) / ( 1 - e_J * cos(E_J*pi/180) );
        }
    }
     xv_J = a_J * (cos(E_J*pi/180) - e_J); //NOT IMPORTANT
     yv_J = a_J * sqrt(1.0 - e_J*e_J) * sin(E_J*pi/180); //NOT IMPORTANT
     v_J = atan2( yv_J, xv_J ); //Mercury's true anomaly
     r_J = sqrt( xv_J*xv_J + yv_J*yv_J ); //Mercury's distance
     xh_J = r_J * ( cos(N_J*pi/180) * cos(v_J+w_J*pi/180) - sin(N_J*pi/180) * sin(v_J+w_J*pi/180) * cos(i_J*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_J = r_J * ( sin(N_J*pi/180) * cos(v_J+w_J*pi/180) + cos(N_J*pi/180) * sin(v_J+w_J*pi/180) * cos(i_J*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_J = r_J * ( sin(v_J+w_J*pi/180) * sin(i_J*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_J = atan2( yh_J, xh_J )*180/pi; //Mercury's ecliptic longitude
     latecl_J = atan2( zh_J, sqrt(xh_J*xh_J+yh_J*yh_J) )*180/pi; //Mercury's ecliptic latitude
	 M_SA=316.9670 + 0.0334442282 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_S>(360))
		M_S = M_S - 360;
	while(M_S<0)
		M_S = M_S + 360;
	 M_U=142.5905 + 0.011725806 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_U>(360))
		M_U = M_U - 360;
	while(M_U<0)
		M_U = M_U + 360;
	
    
	 lonecl_J = lonecl_J -0.332 * sin(2*M_J*pi/180 - 5*M_SA*pi/180 - 67.6*pi/180 )
    -0.056 * sin(2*M_J*pi/180 - 2*M_SA*pi/180 + 21*pi/180 )
    +0.042 * sin(3*M_J*pi/180 - 5*M_SA*pi/180 + 21*pi/180 )
    -0.036 * sin(M_J*pi/180 - 2*M_SA*pi/180)
    +0.022 * cos(M_J*pi/180 - M_SA*pi/180)
    +0.023 * sin(2*M_J*pi/180 - 3*M_SA + 52*pi/180 )
    -0.016 * sin(M_J*pi/180 - 5*M_SA*pi/180 - 69*pi/180 );

	lonecl_SA = lonecl_SA +0.812 * sin(2*M_J*pi/180 - 5*M_SA*pi/180 - 67.6*pi/180 )
    -0.229 * cos(2*M_J*pi/180 - 4*M_SA*pi/180 - 2*pi/180 )
    +0.119 * sin(M_J*pi/180 - 2*M_SA*pi/180 - 3*pi/180 )
    +0.046 * sin(2*M_J*pi/180 - 6*M_SA*pi/180 - 69*pi/180 )
    +0.014 * sin(M_J*pi/180 - 3*M_SA*pi/180 + 32*pi/180 );

	latecl_SA = latecl_SA -0.020 * cos(2*M_J*pi/180 - 4*M_SA*pi/180 - 2*pi/180 )
    +0.018 * sin(2*M_J*pi/180 - 6*M_SA*pi/180 - 49*pi/180 );

	lonecl_U = lonecl_U +0.040 * sin(M_SA*pi/180 - 2*M_U*pi/180 + 6*pi/180 )
    +0.035 * sin(M_SA*pi/180 - 3*M_U*pi/180 + 33*pi/180 )
    -0.015 * sin(M_J*pi/180 - M_U*pi/180 + 20*pi/180 );

	 xg_J = r_J * cos(lonecl_J*pi/180) * cos(latecl_J*pi/180) +xs; //geocentric position
	 yg_J = r_J * sin(lonecl_J*pi/180) * cos(latecl_J*pi/180) +ys; //geocentric position
	 zg_J = r_J * sin(latecl_J*pi/180); //geocentric position
	 xe_J = xg_J; //equatorial coordinates
     ye_J = yg_J * cos(ecl*pi/180) - zg_J * sin(ecl*pi/180); //equatorial coordinates
     ze_J = yg_J * sin(ecl*pi/180) + zg_J * cos(ecl*pi/180); //equatorial coordinates
     RA_J  = atan2( ye_J, xe_J )*180/pi; //Right Ascension
     Dec_J = atan2( ze_J, sqrt(xe_J*xe_J+ye_J*ye_J) )*180/pi; //Declination
	 par_J = (8.794/3600) / r_J;
	 HA_J = LST - RA_J; //hour angle of the Moon
	 x_J = cos(HA_J*pi/180) * cos(Dec_J*pi/180);
     y_J = sin(HA_J*pi/180) * cos(Dec_J*pi/180);
     z_J = sin(Dec_J*pi/180);
	 xhor_J = x_J * sin(35*pi/180) - z_J * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_J = y_J; //Horizontal coordinates
     zhor_J = x_J * cos(35*pi/180) + z_J * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_J  = atan2( yhor_J, xhor_J )*180/pi + 180; //local azimuth
     alt_J = atan2( zhor_J, sqrt(xhor_J*xhor_J+yhor_J*yhor_J) )*180/pi; //local altitude
	 alt_topoc_J = alt_J - par_J * cos(alt_J*pi/180); //topocentric altitude
	 gclat_J = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_J   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_J = atan( tan(gclat_J*pi/180) / cos(HA_J*pi/180) ); //an auxiliary angle
	 topRA_J   = RA_J  - par_J * rho_J * cos(gclat_J*pi/180) * sin(HA_J*pi/180) / cos(Dec_J*pi/180);
     topDecl_J = Dec_J - par_J * rho_J * sin(gclat_J*pi/180) * sin(g_J - Dec_J*pi/180) / sin(g_J);
	 
	  UT_J_in_south = ( topRA_J - GMST0 - 51.5 );
	while(UT_J_in_south>360)
		UT_J_in_south = UT_J_in_south - 360;
	while(UT_J_in_south<0)
		UT_J_in_south = UT_J_in_south +360;
	 cosLHA_J = (sin((-0.583)*pi/180) - sin(35.0*pi/180) * sin(topDecl_J*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_J*pi/180));

	if(cosLHA_J<-1.0){
		cin>>D;
	}
	if(cosLHA_J>1.0){
		cin>>D;
	}
	 set_J;
	rise_J;
	if(cosLHA_J<1 && cosLHA_J>-1)
	{
		LHA_J = acos(cosLHA_J)*180/pi; //Time in degrees
		if(LHA_J>(180))
			LHA_J = LHA_J- 360;
		if(LHA_J<-180)
			LHA_J = LHA_J + 360;
		set_J = UT_J_in_south + LHA_J;
	}
	while(set_J>(360))
		set_J = set_J - 360;
	while(set_J<0)
		set_J = set_J + 360;
	set_J = set_J/15.04107  +4.5;
	perturbed_set_time_J = set_J/24;
	}
	perturbed_set_time_J = perturbed_set_time_J*24;

	
	//////////////////////////////////////////////////////////////////////////////////PERTURBING RISE

	for(int rise=0;rise<4;rise++){
	d = reald ; //day
	d = d + perturbed_rise_time_M; //perturbing time
	N_J=100.4542 + 2.76854E-5 * d;  //longitude of the ascending node
    i_J=1.3030 - 1.557E-7 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    w_J=273.8777 + 1.64505E-5 * d;  //argument of perihelion
    a_J=5.20256  ;  //semi-major axis, or mean distance from Sun
    e_J=0.048498 + 4.469E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M_J=19.8950 + 0.0830853001 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_J>(360))
		M_J = M_J - 360;
	while(M_J<0)
		M_J = M_J + 360;
    
	  E_J = M_J + e_J*(180/pi) * sin(M_J*pi/180) * ( 1.0 + e_J * cos(M_J*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_J>0.05){
        double E0_J=E_J;
        double E1_J=E_J;
        while(fabs(E0_J-E1_J)>0.001){
            E0_J=E1_J;
            E1_J = E_J - ( E_J - e_J*(180/pi) * sin(E_J*pi/180) - M_J ) / ( 1 - e_J * cos(E_J*pi/180) );
        }
    }
     xv_J = a_J * (cos(E_J*pi/180) - e_J); //NOT IMPORTANT
     yv_J = a_J * sqrt(1.0 - e_J*e_J) * sin(E_J*pi/180); //NOT IMPORTANT
     v_J = atan2( yv_J, xv_J ); //Mercury's true anomaly
     r_J = sqrt( xv_J*xv_J + yv_J*yv_J ); //Mercury's distance
     xh_J = r_J * ( cos(N_J*pi/180) * cos(v_J+w_J*pi/180) - sin(N_J*pi/180) * sin(v_J+w_J*pi/180) * cos(i_J*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_J = r_J * ( sin(N_J*pi/180) * cos(v_J+w_J*pi/180) + cos(N_J*pi/180) * sin(v_J+w_J*pi/180) * cos(i_J*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_J = r_J * ( sin(v_J+w_J*pi/180) * sin(i_J*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_J = atan2( yh_J, xh_J )*180/pi; //Mercury's ecliptic longitude
     latecl_J = atan2( zh_J, sqrt(xh_J*xh_J+yh_J*yh_J) )*180/pi; //Mercury's ecliptic latitude
	 M_SA=316.9670 + 0.0334442282 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_S>(360))
		M_S = M_S - 360;
	while(M_S<0)
		M_S = M_S + 360;
	 M_U=142.5905 + 0.011725806 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_U>(360))
		M_U = M_U - 360;
	while(M_U<0)
		M_U = M_U + 360;
	
    
	 lonecl_J = lonecl_J -0.332 * sin(2*M_J*pi/180 - 5*M_SA*pi/180 - 67.6*pi/180 )
    -0.056 * sin(2*M_J*pi/180 - 2*M_SA*pi/180 + 21*pi/180 )
    +0.042 * sin(3*M_J*pi/180 - 5*M_SA*pi/180 + 21*pi/180 )
    -0.036 * sin(M_J*pi/180 - 2*M_SA*pi/180)
    +0.022 * cos(M_J*pi/180 - M_SA*pi/180)
    +0.023 * sin(2*M_J*pi/180 - 3*M_SA + 52*pi/180 )
    -0.016 * sin(M_J*pi/180 - 5*M_SA*pi/180 - 69*pi/180 );

	lonecl_SA = lonecl_SA +0.812 * sin(2*M_J*pi/180 - 5*M_SA*pi/180 - 67.6*pi/180 )
    -0.229 * cos(2*M_J*pi/180 - 4*M_SA*pi/180 - 2*pi/180 )
    +0.119 * sin(M_J*pi/180 - 2*M_SA*pi/180 - 3*pi/180 )
    +0.046 * sin(2*M_J*pi/180 - 6*M_SA*pi/180 - 69*pi/180 )
    +0.014 * sin(M_J*pi/180 - 3*M_SA*pi/180 + 32*pi/180 );

	latecl_SA = latecl_SA -0.020 * cos(2*M_J*pi/180 - 4*M_SA*pi/180 - 2*pi/180 )
    +0.018 * sin(2*M_J*pi/180 - 6*M_SA*pi/180 - 49*pi/180 );

	lonecl_U = lonecl_U +0.040 * sin(M_SA*pi/180 - 2*M_U*pi/180 + 6*pi/180 )
    +0.035 * sin(M_SA*pi/180 - 3*M_U*pi/180 + 33*pi/180 )
    -0.015 * sin(M_J*pi/180 - M_U*pi/180 + 20*pi/180 );

	 xg_J = r_J * cos(lonecl_J*pi/180) * cos(latecl_J*pi/180) +xs; //geocentric position
	 yg_J = r_J * sin(lonecl_J*pi/180) * cos(latecl_J*pi/180) +ys; //geocentric position
	 zg_J = r_J * sin(latecl_J*pi/180); //geocentric position
	 xe_J = xg_J; //equatorial coordinates
     ye_J = yg_J * cos(ecl*pi/180) - zg_J * sin(ecl*pi/180); //equatorial coordinates
     ze_J = yg_J * sin(ecl*pi/180) + zg_J * cos(ecl*pi/180); //equatorial coordinates
     RA_J  = atan2( ye_J, xe_J )*180/pi; //Right Ascension
     Dec_J = atan2( ze_J, sqrt(xe_J*xe_J+ye_J*ye_J) )*180/pi; //Declination
	 par_J = (8.794/3600) / r_J;
	 HA_J = LST - RA_J; //hour angle of the Moon
	 x_J = cos(HA_J*pi/180) * cos(Dec_J*pi/180);
     y_J = sin(HA_J*pi/180) * cos(Dec_J*pi/180);
     z_J = sin(Dec_J*pi/180);
	 xhor_J = x_J * sin(35*pi/180) - z_J * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_J = y_J; //Horizontal coordinates
     zhor_J = x_J * cos(35*pi/180) + z_J * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_J  = atan2( yhor_J, xhor_J )*180/pi + 180; //local azimuth
     alt_J = atan2( zhor_J, sqrt(xhor_J*xhor_J+yhor_J*yhor_J) )*180/pi; //local altitude
	 alt_topoc_J = alt_J - par_J * cos(alt_J*pi/180); //topocentric altitude
	 gclat_J = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_J   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_J = atan( tan(gclat_J*pi/180) / cos(HA_J*pi/180) ); //an auxiliary angle
	 topRA_J   = RA_J  - par_J * rho_J * cos(gclat_J*pi/180) * sin(HA_J*pi/180) / cos(Dec_J*pi/180);
     topDecl_J = Dec_J - par_J * rho_J * sin(gclat_J*pi/180) * sin(g_J - Dec_J*pi/180) / sin(g_J);
	
	  UT_J_in_south = ( topRA_J - GMST0 - 51.5 );
	while(UT_J_in_south>360)
		UT_J_in_south = UT_J_in_south - 360;
	while(UT_J_in_south<0)
		UT_J_in_south = UT_J_in_south +360;
	 cosLHA_J = (sin(-0.583*pi/180) - sin(35.0*pi/180) * sin(topDecl_J*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_J*pi/180));

	if(cosLHA_J<-1.0){
		cin>>D;
	}
	if(cosLHA_J>1.0){
		cin>>D;
	}
	if(cosLHA_J<1 && cosLHA_J>-1)
	{
		LHA_J = acos(cosLHA_J)*180/pi; //Time in degrees
		if(LHA_J>(180))
			LHA_J = LHA_J- 360;
		if(LHA_J<-180)
			LHA_J = LHA_J + 360;
		rise_J = UT_J_in_south - LHA_J;
	}
	while(rise_J>(360))
		rise_J = rise_J - 360;
	while(rise_J<0)
		rise_J = rise_J + 360;
	rise_J = rise_J/15.04107  + 4.5;
	perturbed_rise_time_J = rise_J/24;
	}
	perturbed_rise_time_J = perturbed_rise_time_J*24;
	UT_J_in_south = UT_J_in_south/15;
	UT_J_in_south = UT_J_in_south + 4.5;
	double meridianhour_J = (int)UT_J_in_south;
	double meridianminute_J = (int)((UT_J_in_south - meridianhour_J) * 60);
	double meridiansecond_J = (int)((((UT_J_in_south - meridianhour_J) * 60) - meridianminute_J) * 60);
	double sethour_J = (int)set_J;
	double setminute_J = (int)((set_J - sethour_J) * 60);
	double setsecond_J = (int)((((set_J - sethour_J) * 60) - setminute_J) * 60);
	double risehour_J = (int)rise_J;
	double riseminute_J = (int)((rise_J - risehour_J) * 60);
	double risesecond_J = (int)((((rise_J - risehour_J) * 60) - riseminute_J) * 60);
	cout<<"Jupiter:"<<endl<< "setting time: "<<sethour_J<<" : "<<setminute_J<<" : "<<setsecond_J<<endl<<"rising time: "<<risehour_J<<" : "<<riseminute_J<<" : "<<risesecond_J<<endl<<"meridian time:"<<meridianhour_J<<" : "<<meridianminute_J<<" : "<<meridiansecond_J<<endl;

	
	////////////////////////////////Saturn's rise and set times////////////////////////////////
      	d = reald + 0.3125;
	N_SA=113.6634 + 2.38980E-5 * d;  //longitude of the ascending node
    i_SA=2.4886 - 1.081E-7 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    w_SA=339.3939 + 2.97661E-5 * d;  //argument of perihelion
    a_SA=9.55475  ;  //semi-major axis, or mean distance from Sun
    e_SA=0.055546 - 9.499E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M_SA=316.9670 + 0.0334442282 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_S>(360))
		M_S = M_S - 360;
	while(M_S<0)
		M_S = M_S + 360;
    
	  E_SA = M_SA + e_SA*(180/pi) * sin(M_SA*pi/180) * ( 1.0 + e_SA * cos(M_SA*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_SA>0.05){
        double E0_SA=E_SA;
        double E1_SA=E_SA;
        while(fabs(E0_SA-E1_SA)>0.001){
            E0_SA=E1_SA;
            E1_SA = E_SA - ( E_SA - e_SA*(180/pi) * sin(E_SA*pi/180) - M_SA ) / ( 1 - e_SA * cos(E_SA*pi/180) );
        }
    }
     xv_SA = a_SA * (cos(E_SA*pi/180) - e_SA); //NOT IMPORTANT
     yv_SA = a_SA * sqrt(1.0 - e_SA*e_SA) * sin(E_SA*pi/180); //NOT IMPORTANT
     v_SA = atan2( yv_SA, xv_SA ); //Mercury's true anomaly
     r_SA = sqrt( xv_SA*xv_SA + yv_SA*yv_SA ); //Mercury's distance
     xh_SA = r_SA * ( cos(N_SA*pi/180) * cos(v_SA+w_SA*pi/180) - sin(N_SA*pi/180) * sin(v_SA+w_SA*pi/180) * cos(i_SA*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_SA = r_SA * ( sin(N_SA*pi/180) * cos(v_SA+w_SA*pi/180) + cos(N_SA*pi/180) * sin(v_SA+w_SA*pi/180) * cos(i_SA*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_SA = r_SA * ( sin(v_SA+w_SA*pi/180) * sin(i_SA*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_SA = atan2( yh_SA, xh_SA )*180/pi; //Mercury's ecliptic longitude
     latecl_SA = atan2( zh_SA, sqrt(xh_SA*xh_SA+yh_SA*yh_SA) )*180/pi; //Mercury's ecliptic latitude
	 M_J=19.8950 + 0.0830853001 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_J>(360))
		M_J = M_J - 360;
	while(M_J<0)
		M_J = M_J + 360;
     M_U=142.5905 + 0.011725806 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_U>(360))
		M_U = M_U - 360;
	while(M_U<0)
		M_U = M_U + 360;
	
	 lonecl_J = lonecl_J -0.332 * sin(2*M_J*pi/180 - 5*M_SA*pi/180 - 67.6*pi/180 )
    -0.056 * sin(2*M_J*pi/180 - 2*M_SA*pi/180 + 21*pi/180 )
    +0.042 * sin(3*M_J*pi/180 - 5*M_SA*pi/180 + 21*pi/180 )
    -0.036 * sin(M_J*pi/180 - 2*M_SA*pi/180)
    +0.022 * cos(M_J*pi/180 - M_SA*pi/180)
    +0.023 * sin(2*M_J*pi/180 - 3*M_SA + 52*pi/180 )
    -0.016 * sin(M_J*pi/180 - 5*M_SA*pi/180 - 69*pi/180 );

	lonecl_SA = lonecl_SA +0.812 * sin(2*M_J*pi/180 - 5*M_SA*pi/180 - 67.6*pi/180 )
    -0.229 * cos(2*M_J*pi/180 - 4*M_SA*pi/180 - 2*pi/180 )
    +0.119 * sin(M_J*pi/180 - 2*M_SA*pi/180 - 3*pi/180 )
    +0.046 * sin(2*M_J*pi/180 - 6*M_SA*pi/180 - 69*pi/180 )
    +0.014 * sin(M_J*pi/180 - 3*M_SA*pi/180 + 32*pi/180 );

	latecl_SA = latecl_SA -0.020 * cos(2*M_J*pi/180 - 4*M_SA*pi/180 - 2*pi/180 )
    +0.018 * sin(2*M_J*pi/180 - 6*M_SA*pi/180 - 49*pi/180 );

	lonecl_U = lonecl_U +0.040 * sin(M_SA*pi/180 - 2*M_U*pi/180 + 6*pi/180 )
    +0.035 * sin(M_SA*pi/180 - 3*M_U*pi/180 + 33*pi/180 )
    -0.015 * sin(M_J*pi/180 - M_U*pi/180 + 20*pi/180 );
	 xg_SA = r_SA * cos(lonecl_SA*pi/180) * cos(latecl_SA*pi/180) +xs; //geocentric position
	 yg_SA = r_SA * sin(lonecl_SA*pi/180) * cos(latecl_SA*pi/180) +ys; //geocentric position
	 zg_SA = r_SA * sin(latecl_SA*pi/180); //geocentric position
	 xe_SA = xg_SA; //equatorial coordinates
     ye_SA = yg_SA * cos(ecl*pi/180) - zg_SA * sin(ecl*pi/180); //equatorial coordinates
     ze_SA = yg_SA * sin(ecl*pi/180) + zg_SA * cos(ecl*pi/180); //equatorial coordinates
     RA_SA  = atan2( ye_SA, xe_SA )*180/pi; //Right Ascension
     Dec_SA = atan2( ze_SA, sqrt(xe_SA*xe_SA+ye_SA*ye_SA) )*180/pi; //Declination
	 par_SA = (8.794/3600) / r_SA;
	 HA_SA = LST - RA_SA; //hour angle of the Moon
	 x_SA = cos(HA_SA*pi/180) * cos(Dec_SA*pi/180);
     y_SA = sin(HA_SA*pi/180) * cos(Dec_SA*pi/180);
     z_SA = sin(Dec_SA*pi/180);
	 xhor_SA = x_SA * sin(35*pi/180) - z_SA * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_SA = y_SA; //Horizontal coordinates
     zhor_SA = x_SA * cos(35*pi/180) + z_SA * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_SA  = atan2( yhor_SA, xhor_SA )*180/pi + 180; //local azimuth
     alt_SA = atan2( zhor_SA, sqrt(xhor_SA*xhor_SA+yhor_SA*yhor_SA) )*180/pi; //local altitude
	 alt_topoc_SA = alt_SA - par_SA * cos(alt_SA*pi/180); //topocentric altitude
	 gclat_SA = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_SA   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_SA = atan( tan(gclat_SA*pi/180) / cos(HA_SA*pi/180) ); //an auxiliary angle
	 topRA_SA   = RA_SA  - par_SA * rho_SA * cos(gclat_SA*pi/180) * sin(HA_SA*pi/180) / cos(Dec_SA*pi/180);
     topDecl_SA = Dec_SA - par_SA * rho_SA * sin(gclat_SA*pi/180) * sin(g_SA - Dec_SA*pi/180) / sin(g_SA);
	
	 double UT_SA_in_south = ( topRA_SA - GMST0 - 51.5 );
	while(UT_SA_in_south>360)
		UT_SA_in_south = UT_SA_in_south - 360;
	while(UT_SA_in_south<0)
		UT_SA_in_south = UT_SA_in_south +360;
	double cosLHA_SA = (sin((-0.583)*pi/180) - sin(35.0*pi/180) * sin(topDecl_SA*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_SA*pi/180));
	double LHA_SA;
	if(cosLHA_SA<-1.0){
		cin>>D;
	}
	if(cosLHA_SA>1.0){
		cin>>D;
	}
	double set_SA;
	double rise_SA;
	if(cosLHA_SA<1 && cosLHA_SA>-1)
	{
		LHA_SA = acos(cosLHA_SA)*180/pi; //Time in degrees
		if(LHA_SA>(180))
			LHA_SA = LHA_SA- 360;
		if(LHA_SA<-180)
			LHA_SA = LHA_SA + 360;
		set_SA = UT_SA_in_south + LHA_SA;
		rise_SA = UT_SA_in_south - LHA_SA;
	}
	while(set_SA>(360))
		set_SA = set_SA - 360;
	while(set_SA<0)
		set_SA = set_SA + 360;
	while(rise_SA>(360))
		rise_SA = rise_SA - 360;
	while(rise_SA<0)
		rise_SA = rise_SA + 360;
	set_SA = set_SA/15.04107  +4.5;
	rise_SA = rise_SA/15.04107  +4.5;
	//////////////////////////////////////////////////////////////////////////////////PERTURBING SET
	double perturbed_rise_time_SA = rise_SA/24;
	double perturbed_set_time_SA = set_SA/24;

	for(int set=0;set<4;set++){
	d = reald ; //day
	d = d + perturbed_set_time_SA; //perturbing time
	N_SA=113.6634 + 2.38980E-5 * d;  //longitude of the ascending node
    i_SA=2.4886 - 1.081E-7 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    w_SA=339.3939 + 2.97661E-5 * d;  //argument of perihelion
    a_SA=9.55475  ;  //semi-major axis, or mean distance from Sun
    e_SA=0.055546 - 9.499E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M_SA=316.9670 + 0.0334442282 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_S>(360))
		M_S = M_S - 360;
	while(M_S<0)
		M_S = M_S + 360;
    
	  E_SA = M_SA + e_SA*(180/pi) * sin(M_SA*pi/180) * ( 1.0 + e_SA * cos(M_SA*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_SA>0.05){
        double E0_SA=E_SA;
        double E1_SA=E_SA;
        while(fabs(E0_SA-E1_SA)>0.001){
            E0_SA=E1_SA;
            E1_SA = E_SA - ( E_SA - e_SA*(180/pi) * sin(E_SA*pi/180) - M_SA ) / ( 1 - e_SA * cos(E_SA*pi/180) );
        }
    }
     xv_SA = a_SA * (cos(E_SA*pi/180) - e_SA); //NOT IMPORTANT
     yv_SA = a_SA * sqrt(1.0 - e_SA*e_SA) * sin(E_SA*pi/180); //NOT IMPORTANT
     v_SA = atan2( yv_SA, xv_SA ); //Mercury's true anomaly
     r_SA = sqrt( xv_SA*xv_SA + yv_SA*yv_SA ); //Mercury's distance
     xh_SA = r_SA * ( cos(N_SA*pi/180) * cos(v_SA+w_SA*pi/180) - sin(N_SA*pi/180) * sin(v_SA+w_SA*pi/180) * cos(i_SA*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_SA = r_SA * ( sin(N_SA*pi/180) * cos(v_SA+w_SA*pi/180) + cos(N_SA*pi/180) * sin(v_SA+w_SA*pi/180) * cos(i_SA*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_SA = r_SA * ( sin(v_SA+w_SA*pi/180) * sin(i_SA*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_SA = atan2( yh_SA, xh_SA )*180/pi; //Mercury's ecliptic longitude
     latecl_SA = atan2( zh_SA, sqrt(xh_SA*xh_SA+yh_SA*yh_SA) )*180/pi; //Mercury's ecliptic latitude
	M_J=19.8950 + 0.0830853001 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_J>(360))
		M_J = M_J - 360;
	while(M_J<0)
		M_J = M_J + 360;
	 M_U=142.5905 + 0.011725806 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_U>(360))
		M_U = M_U - 360;
	while(M_U<0)
		M_U = M_U + 360;
	
    
	lonecl_J = lonecl_J -0.332 * sin(2*M_J*pi/180 - 5*M_SA*pi/180 - 67.6*pi/180 )
    -0.056 * sin(2*M_J*pi/180 - 2*M_SA*pi/180 + 21*pi/180 )
    +0.042 * sin(3*M_J*pi/180 - 5*M_SA*pi/180 + 21*pi/180 )
    -0.036 * sin(M_J*pi/180 - 2*M_SA*pi/180)
    +0.022 * cos(M_J*pi/180 - M_SA*pi/180)
    +0.023 * sin(2*M_J*pi/180 - 3*M_SA + 52*pi/180 )
    -0.016 * sin(M_J*pi/180 - 5*M_SA*pi/180 - 69*pi/180 );

	lonecl_SA = lonecl_SA +0.812 * sin(2*M_J*pi/180 - 5*M_SA*pi/180 - 67.6*pi/180 )
    -0.229 * cos(2*M_J*pi/180 - 4*M_SA*pi/180 - 2*pi/180 )
    +0.119 * sin(M_J*pi/180 - 2*M_SA*pi/180 - 3*pi/180 )
    +0.046 * sin(2*M_J*pi/180 - 6*M_SA*pi/180 - 69*pi/180 )
    +0.014 * sin(M_J*pi/180 - 3*M_SA*pi/180 + 32*pi/180 );

	latecl_SA = latecl_SA -0.020 * cos(2*M_J*pi/180 - 4*M_SA*pi/180 - 2*pi/180 )
    +0.018 * sin(2*M_J*pi/180 - 6*M_SA*pi/180 - 49*pi/180 );

	lonecl_U = lonecl_U +0.040 * sin(M_SA*pi/180 - 2*M_U*pi/180 + 6*pi/180 )
    +0.035 * sin(M_SA*pi/180 - 3*M_U*pi/180 + 33*pi/180 )
    -0.015 * sin(M_J*pi/180 - M_U*pi/180 + 20*pi/180 );
	 xg_SA = r_SA * cos(lonecl_SA*pi/180) * cos(latecl_SA*pi/180) +xs; //geocentric position
	 yg_SA = r_SA * sin(lonecl_SA*pi/180) * cos(latecl_SA*pi/180) +ys; //geocentric position
	 zg_SA = r_SA * sin(latecl_SA*pi/180); //geocentric position
	 xe_SA = xg_SA; //equatorial coordinates
     ye_SA = yg_SA * cos(ecl*pi/180) - zg_SA * sin(ecl*pi/180); //equatorial coordinates
     ze_SA = yg_SA * sin(ecl*pi/180) + zg_SA * cos(ecl*pi/180); //equatorial coordinates
     RA_SA  = atan2( ye_SA, xe_SA )*180/pi; //Right Ascension
     Dec_SA = atan2( ze_SA, sqrt(xe_SA*xe_SA+ye_SA*ye_SA) )*180/pi; //Declination
	 par_SA = (8.794/3600) / r_SA;
	 HA_SA = LST - RA_SA; //hour angle of the Moon
	 x_SA = cos(HA_SA*pi/180) * cos(Dec_SA*pi/180);
     y_SA = sin(HA_SA*pi/180) * cos(Dec_SA*pi/180);
     z_SA = sin(Dec_SA*pi/180);
	 xhor_SA = x_SA * sin(35*pi/180) - z_SA * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_SA = y_SA; //Horizontal coordinates
     zhor_SA = x_SA * cos(35*pi/180) + z_SA * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_SA  = atan2( yhor_SA, xhor_SA )*180/pi + 180; //local azimuth
     alt_SA = atan2( zhor_SA, sqrt(xhor_SA*xhor_SA+yhor_SA*yhor_SA) )*180/pi; //local altitude
	 alt_topoc_SA = alt_SA - par_SA * cos(alt_SA*pi/180); //topocentric altitude
	 gclat_SA = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_SA   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_SA = atan( tan(gclat_SA*pi/180) / cos(HA_SA*pi/180) ); //an auxiliary angle
	 topRA_SA   = RA_SA  - par_SA * rho_SA * cos(gclat_SA*pi/180) * sin(HA_SA*pi/180) / cos(Dec_SA*pi/180);
     topDecl_SA = Dec_SA - par_SA * rho_SA * sin(gclat_SA*pi/180) * sin(g_SA - Dec_SA*pi/180) / sin(g_SA);
	 
	  UT_SA_in_south = ( topRA_SA - GMST0 - 51.5 );
	while(UT_SA_in_south>360)
		UT_SA_in_south = UT_SA_in_south - 360;
	while(UT_SA_in_south<0)
		UT_SA_in_south = UT_SA_in_south +360;
	 cosLHA_SA = (sin((-0.583)*pi/180) - sin(35.0*pi/180) * sin(topDecl_SA*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_SA*pi/180));

	if(cosLHA_SA<-1.0){
		cin>>D;
	}
	if(cosLHA_SA>1.0){
		cin>>D;
	}
	 set_SA;
	rise_SA;
	if(cosLHA_SA<1 && cosLHA_SA>-1)
	{
		LHA_SA = acos(cosLHA_SA)*180/pi; //Time in degrees
		if(LHA_SA>(180))
			LHA_SA = LHA_SA- 360;
		if(LHA_SA<-180)
			LHA_SA = LHA_SA + 360;
		set_SA = UT_SA_in_south + LHA_SA;
	}
	while(set_SA>(360))
		set_SA = set_SA - 360;
	while(set_SA<0)
		set_SA = set_SA + 360;
	set_SA = set_SA/15.04107  +4.5;
	perturbed_set_time_SA = set_SA/24;
	}
	perturbed_set_time_SA = perturbed_set_time_SA*24;

	
	//////////////////////////////////////////////////////////////////////////////////PERTURBING RISE

	for(int rise=0;rise<4;rise++){
	d = reald ; //day
	d = d + perturbed_rise_time_M; //perturbing time
	N_SA=113.6634 + 2.38980E-5 * d;  //longitude of the ascending node
    i_SA=2.4886 - 1.081E-7 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    w_SA=339.3939 + 2.97661E-5 * d;  //argument of perihelion
    a_SA=9.55475  ;  //semi-major axis, or mean distance from Sun
    e_SA=0.055546 - 9.499E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M_SA=316.9670 + 0.0334442282 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_S>(360))
		M_S = M_S - 360;
	while(M_S<0)
		M_S = M_S + 360;
    
	  E_SA = M_SA + e_SA*(180/pi) * sin(M_SA*pi/180) * ( 1.0 + e_SA * cos(M_SA*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_SA>0.05){
        double E0_SA=E_SA;
        double E1_SA=E_SA;
        while(fabs(E0_SA-E1_SA)>0.001){
            E0_SA=E1_SA;
            E1_SA = E_SA - ( E_SA - e_SA*(180/pi) * sin(E_SA*pi/180) - M_SA ) / ( 1 - e_SA * cos(E_SA*pi/180) );
        }
    }
     xv_SA = a_SA * (cos(E_SA*pi/180) - e_SA); //NOT IMPORTANT
     yv_SA = a_SA * sqrt(1.0 - e_SA*e_SA) * sin(E_SA*pi/180); //NOT IMPORTANT
     v_SA = atan2( yv_SA, xv_SA ); //Mercury's true anomaly
     r_SA = sqrt( xv_SA*xv_SA + yv_SA*yv_SA ); //Mercury's distance
     xh_SA = r_SA * ( cos(N_SA*pi/180) * cos(v_SA+w_SA*pi/180) - sin(N_SA*pi/180) * sin(v_SA+w_SA*pi/180) * cos(i_SA*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_SA = r_SA * ( sin(N_SA*pi/180) * cos(v_SA+w_SA*pi/180) + cos(N_SA*pi/180) * sin(v_SA+w_SA*pi/180) * cos(i_SA*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_SA = r_SA * ( sin(v_SA+w_SA*pi/180) * sin(i_SA*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_SA = atan2( yh_SA, xh_SA )*180/pi; //Mercury's ecliptic longitude
     latecl_SA = atan2( zh_SA, sqrt(xh_SA*xh_SA+yh_SA*yh_SA) )*180/pi; //Mercury's ecliptic latitude
	 M_J=19.8950 + 0.0830853001 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_J>(360))
		M_J = M_J - 360;
	while(M_J<0)
		M_J = M_J + 360;
	 M_U=142.5905 + 0.011725806 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_U>(360))
		M_U = M_U - 360;
	while(M_U<0)
		M_U = M_U + 360;
	
    
	 lonecl_J = lonecl_J -0.332 * sin(2*M_J*pi/180 - 5*M_SA*pi/180 - 67.6*pi/180 )
    -0.056 * sin(2*M_J*pi/180 - 2*M_SA*pi/180 + 21*pi/180 )
    +0.042 * sin(3*M_J*pi/180 - 5*M_SA*pi/180 + 21*pi/180 )
    -0.036 * sin(M_J*pi/180 - 2*M_SA*pi/180)
    +0.022 * cos(M_J*pi/180 - M_SA*pi/180)
    +0.023 * sin(2*M_J*pi/180 - 3*M_SA + 52*pi/180 )
    -0.016 * sin(M_J*pi/180 - 5*M_SA*pi/180 - 69*pi/180 );

	lonecl_SA = lonecl_SA +0.812 * sin(2*M_J*pi/180 - 5*M_SA*pi/180 - 67.6*pi/180 )
    -0.229 * cos(2*M_J*pi/180 - 4*M_SA*pi/180 - 2*pi/180 )
    +0.119 * sin(M_J*pi/180 - 2*M_SA*pi/180 - 3*pi/180 )
    +0.046 * sin(2*M_J*pi/180 - 6*M_SA*pi/180 - 69*pi/180 )
    +0.014 * sin(M_J*pi/180 - 3*M_SA*pi/180 + 32*pi/180 );

	latecl_SA = latecl_SA -0.020 * cos(2*M_J*pi/180 - 4*M_SA*pi/180 - 2*pi/180 )
    +0.018 * sin(2*M_J*pi/180 - 6*M_SA*pi/180 - 49*pi/180 );

	lonecl_U = lonecl_U +0.040 * sin(M_SA*pi/180 - 2*M_U*pi/180 + 6*pi/180 )
    +0.035 * sin(M_SA*pi/180 - 3*M_U*pi/180 + 33*pi/180 )
    -0.015 * sin(M_J*pi/180 - M_U*pi/180 + 20*pi/180 );
	 xg_SA = r_SA * cos(lonecl_SA*pi/180) * cos(latecl_SA*pi/180) +xs; //geocentric position
	 yg_SA = r_SA * sin(lonecl_SA*pi/180) * cos(latecl_SA*pi/180) +ys; //geocentric position
	 zg_SA = r_SA * sin(latecl_SA*pi/180); //geocentric position
	 xe_SA = xg_SA; //equatorial coordinates
     ye_SA = yg_SA * cos(ecl*pi/180) - zg_SA * sin(ecl*pi/180); //equatorial coordinates
     ze_SA = yg_SA * sin(ecl*pi/180) + zg_SA * cos(ecl*pi/180); //equatorial coordinates
     RA_SA  = atan2( ye_SA, xe_SA )*180/pi; //Right Ascension
     Dec_SA = atan2( ze_SA, sqrt(xe_SA*xe_SA+ye_SA*ye_SA) )*180/pi; //Declination
	 par_SA = (8.794/3600) / r_SA;
	 HA_SA = LST - RA_SA; //hour angle of the Moon
	 x_SA = cos(HA_SA*pi/180) * cos(Dec_SA*pi/180);
     y_SA = sin(HA_SA*pi/180) * cos(Dec_SA*pi/180);
     z_SA = sin(Dec_SA*pi/180);
	 xhor_SA = x_SA * sin(35*pi/180) - z_SA * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_SA = y_SA; //Horizontal coordinates
     zhor_SA = x_SA * cos(35*pi/180) + z_SA * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_SA  = atan2( yhor_SA, xhor_SA )*180/pi + 180; //local azimuth
     alt_SA = atan2( zhor_SA, sqrt(xhor_SA*xhor_SA+yhor_SA*yhor_SA) )*180/pi; //local altitude
	 alt_topoc_SA = alt_SA - par_SA * cos(alt_SA*pi/180); //topocentric altitude
	 gclat_SA = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_SA   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_SA = atan( tan(gclat_SA*pi/180) / cos(HA_SA*pi/180) ); //an auxiliary angle
	 topRA_SA   = RA_SA  - par_SA * rho_SA * cos(gclat_SA*pi/180) * sin(HA_SA*pi/180) / cos(Dec_SA*pi/180);
     topDecl_SA = Dec_SA - par_SA * rho_SA * sin(gclat_SA*pi/180) * sin(g_SA - Dec_SA*pi/180) / sin(g_SA);
	
	  UT_SA_in_south = ( topRA_SA - GMST0 - 51.5 );
	while(UT_SA_in_south>360)
		UT_SA_in_south = UT_SA_in_south - 360;
	while(UT_SA_in_south<0)
		UT_SA_in_south = UT_SA_in_south +360;
	 cosLHA_SA = (sin(-0.583*pi/180) - sin(35.0*pi/180) * sin(topDecl_SA*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_SA*pi/180));

	if(cosLHA_SA<-1.0){
		cin>>D;
	}
	if(cosLHA_SA>1.0){
		cin>>D;
	}
	if(cosLHA_SA<1 && cosLHA_SA>-1)
	{
		LHA_SA = acos(cosLHA_SA)*180/pi; //Time in degrees
		if(LHA_SA>(180))
			LHA_SA = LHA_SA- 360;
		if(LHA_SA<-180)
			LHA_SA = LHA_SA + 360;
		rise_SA = UT_SA_in_south - LHA_SA;
	}
	while(rise_SA>(360))
		rise_SA = rise_SA - 360;
	while(rise_SA<0)
		rise_SA = rise_SA + 360;
	rise_SA = rise_SA/15.04107  + 4.5;
	perturbed_rise_time_SA = rise_SA/24;
	}
	perturbed_rise_time_SA = perturbed_rise_time_SA*24;
	UT_SA_in_south = UT_SA_in_south/15;
	UT_SA_in_south = UT_SA_in_south + 4.5;
	double meridianhour_SA = (int)UT_SA_in_south;
	double meridianminute_SA = (int)((UT_SA_in_south - meridianhour_SA) * 60);
	double meridiansecond_SA = (int)((((UT_SA_in_south - meridianhour_SA) * 60) - meridianminute_SA) * 60);
	double sethour_SA = (int)set_SA;
	double setminute_SA = (int)((set_SA - sethour_SA) * 60);
	double setsecond_SA = (int)((((set_SA - sethour_SA) * 60) - setminute_SA) * 60);
	double risehour_SA = (int)rise_SA;
	double riseminute_SA = (int)((rise_SA - risehour_SA) * 60);
	double risesecond_SA = (int)((((rise_SA - risehour_SA) * 60) - riseminute_SA) * 60);
	cout<<"Saturn:"<<endl<< "setting time: "<<sethour_SA<<" : "<<setminute_SA<<" : "<<setsecond_SA<<endl<<"rising time: "<<risehour_SA<<" : "<<riseminute_SA<<" : "<<risesecond_SA<<endl<<"meridian time:"<<meridianhour_SA<<" : "<<meridianminute_SA<<" : "<<meridiansecond_SA<<endl;

	////////////////////////////////Uranus's rise and set times////////////////////////////////
    d = reald + 0.3125;
	N_U=74.0005 + 1.3978E-5 * d;  //longitude of the ascending node
    i_U=0.7733 + 1.9E-8 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    w_U=96.6612 + 3.0565E-5 * d;  //argument of perihelion
    a_U=19.18171 - 1.55E-8 * d;  //semi-major axis, or mean distance from Sun
    e_U=0.047318 + 7.45E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M_U=142.5905 + 0.011725806 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_U>(360))
		M_U = M_U - 360;
	while(M_U<0)
		M_U = M_U + 360;
	  E_U = M_U + e_U*(180/pi) * sin(M_U*pi/180) * ( 1.0 + e_U * cos(M_U*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_U>0.05){
        double E0_U=E_U;
        double E1_U=E_U;
        while(fabs(E0_U-E1_U)>0.001){
            E0_U=E1_U;
            E1_U = E_U - ( E_U - e_U*(180/pi) * sin(E_U*pi/180) - M_U ) / ( 1 - e_U * cos(E_U*pi/180) );
        }
    }
     xv_U = a_U * (cos(E_U*pi/180) - e_U); //NOT IMPORTANT
     yv_U = a_U * sqrt(1.0 - e_U*e_U) * sin(E_U*pi/180); //NOT IMPORTANT
     v_U = atan2( yv_U, xv_U ); //Mercury's true anomaly
     r_U = sqrt( xv_U*xv_U + yv_U*yv_U ); //Mercury's distance
     xh_U = r_U * ( cos(N_U*pi/180) * cos(v_U+w_U*pi/180) - sin(N_U*pi/180) * sin(v_U+w_U*pi/180) * cos(i_U*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_U = r_U * ( sin(N_U*pi/180) * cos(v_U+w_U*pi/180) + cos(N_U*pi/180) * sin(v_U+w_U*pi/180) * cos(i_U*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_U = r_U * ( sin(v_U+w_U*pi/180) * sin(i_U*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_U = atan2( yh_U, xh_U )*180/pi; //Mercury's ecliptic longitude
     latecl_U = atan2( zh_U, sqrt(xh_U*xh_U+yh_U*yh_U) )*180/pi; //Mercury's ecliptic latitude
	 xg_U = r_U * cos(lonecl_U*pi/180) * cos(latecl_U*pi/180) +xs; //geocentric position
	 yg_U = r_U * sin(lonecl_U*pi/180) * cos(latecl_U*pi/180) +ys; //geocentric position
	 zg_U = r_U * sin(latecl_U*pi/180); //geocentric position
	 xe_U = xg_U; //equatorial coordinates
     ye_U = yg_U * cos(ecl*pi/180) - zg_U * sin(ecl*pi/180); //equatorial coordinates
     ze_U = yg_U * sin(ecl*pi/180) + zg_U * cos(ecl*pi/180); //equatorial coordinates
     RA_U  = atan2( ye_U, xe_U )*180/pi; //Right Ascension
     Dec_U = atan2( ze_U, sqrt(xe_U*xe_U+ye_U*ye_U) )*180/pi; //Declination
	 par_U = (8.794/3600) / r_U;
	 HA_U = LST - RA_U; //hour angle of the Moon
	 x_U = cos(HA_U*pi/180) * cos(Dec_U*pi/180);
     y_U = sin(HA_U*pi/180) * cos(Dec_U*pi/180);
     z_U = sin(Dec_U*pi/180);
	 xhor_U = x_U * sin(35*pi/180) - z_U * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_U = y_U; //Horizontal coordinates
     zhor_U = x_U * cos(35*pi/180) + z_U * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_U  = atan2( yhor_U, xhor_U )*180/pi + 180; //local azimuth
     alt_U = atan2( zhor_U, sqrt(xhor_U*xhor_U+yhor_U*yhor_U) )*180/pi; //local altitude
	 alt_topoc_U = alt_U - par_U * cos(alt_U*pi/180); //topocentric altitude
	 gclat_U = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_U   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_U = atan( tan(gclat_U*pi/180) / cos(HA_U*pi/180) ); //an auxiliary angle
	 topRA_U   = RA_U  - par_U * rho_U * cos(gclat_U*pi/180) * sin(HA_U*pi/180) / cos(Dec_U*pi/180);
     topDecl_U = Dec_U - par_U * rho_U * sin(gclat_U*pi/180) * sin(g_U - Dec_U*pi/180) / sin(g_U);
	
	 double UT_U_in_south = ( topRA_U - GMST0 - 51.5 );
	while(UT_U_in_south>360)
		UT_U_in_south = UT_U_in_south - 360;
	while(UT_U_in_south<0)
		UT_U_in_south = UT_U_in_south +360;
	double cosLHA_U = (sin((-0.583)*pi/180) - sin(35.0*pi/180) * sin(topDecl_U*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_U*pi/180));
	double LHA_U;
	if(cosLHA_U<-1.0){
		cin>>D;
	}
	if(cosLHA_U>1.0){
		cin>>D;
	}
	double set_U;
	double rise_U;
	if(cosLHA_U<1 && cosLHA_U>-1)
	{
		LHA_U = acos(cosLHA_U)*180/pi; //Time in degrees
		if(LHA_U>(180))
			LHA_U = LHA_U- 360;
		if(LHA_U<-180)
			LHA_U = LHA_U + 360;
		set_U = UT_U_in_south + LHA_U;
		rise_U = UT_U_in_south - LHA_U;
	}
	while(set_U>(360))
		set_U = set_U - 360;
	while(set_U<0)
		set_U = set_U + 360;
	while(rise_U>(360))
		rise_U = rise_U - 360;
	while(rise_U<0)
		rise_U = rise_U + 360;
	set_U = set_U/15.04107  +4.5;
	rise_U = rise_U/15.04107  +4.5;
	//////////////////////////////////////////////////////////////////////////////////PERTURBING SET
	double perturbed_rise_time_U = rise_U/24;
	double perturbed_set_time_U = set_U/24;

	for(int set=0;set<4;set++){
	d = reald ; //day
	d = d + perturbed_set_time_U; //perturbing time
	N_U=74.0005 + 1.3978E-5 * d;  //longitude of the ascending node
    i_U=0.7733 + 1.9E-8 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    w_U=96.6612 + 3.0565E-5 * d;  //argument of perihelion
    a_U=19.18171 - 1.55E-8 * d;  //semi-major axis, or mean distance from Sun
    e_U=0.047318 + 7.45E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M_U=142.5905 + 0.011725806 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_U>(360))
		M_U = M_U - 360;
	while(M_U<0)
		M_U = M_U + 360;
	
	  E_U = M_U + e_U*(180/pi) * sin(M_U*pi/180) * ( 1.0 + e_U * cos(M_U*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_U>0.05){
        double E0_U=E_U;
        double E1_U=E_U;
        while(fabs(E0_U-E1_U)>0.001){
            E0_U=E1_U;
            E1_U = E_U - ( E_U - e_U*(180/pi) * sin(E_U*pi/180) - M_U ) / ( 1 - e_U * cos(E_U*pi/180) );
        }
    }
     xv_U = a_U * (cos(E_U*pi/180) - e_U); //NOT IMPORTANT
     yv_U = a_U * sqrt(1.0 - e_U*e_U) * sin(E_U*pi/180); //NOT IMPORTANT
     v_U = atan2( yv_U, xv_U ); //Mercury's true anomaly
     r_U = sqrt( xv_U*xv_U + yv_U*yv_U ); //Mercury's distance
     xh_U = r_U * ( cos(N_U*pi/180) * cos(v_U+w_U*pi/180) - sin(N_U*pi/180) * sin(v_U+w_U*pi/180) * cos(i_U*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_U = r_U * ( sin(N_U*pi/180) * cos(v_U+w_U*pi/180) + cos(N_U*pi/180) * sin(v_U+w_U*pi/180) * cos(i_U*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_U = r_U * ( sin(v_U+w_U*pi/180) * sin(i_U*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_U = atan2( yh_U, xh_U )*180/pi; //Mercury's ecliptic longitude
     latecl_U = atan2( zh_U, sqrt(xh_U*xh_U+yh_U*yh_U) )*180/pi; //Mercury's ecliptic latitude
	 xg_U = r_U * cos(lonecl_U*pi/180) * cos(latecl_U*pi/180) +xs; //geocentric position
	 yg_U = r_U * sin(lonecl_U*pi/180) * cos(latecl_U*pi/180) +ys; //geocentric position
	 zg_U = r_U * sin(latecl_U*pi/180); //geocentric position
	 xe_U = xg_U; //equatorial coordinates
     ye_U = yg_U * cos(ecl*pi/180) - zg_U * sin(ecl*pi/180); //equatorial coordinates
     ze_U = yg_U * sin(ecl*pi/180) + zg_U * cos(ecl*pi/180); //equatorial coordinates
     RA_U  = atan2( ye_U, xe_U )*180/pi; //Right Ascension
     Dec_U = atan2( ze_U, sqrt(xe_U*xe_U+ye_U*ye_U) )*180/pi; //Declination
	 par_U = (8.794/3600) / r_U;
	 HA_U = LST - RA_U; //hour angle of the Moon
	 x_U = cos(HA_U*pi/180) * cos(Dec_U*pi/180);
     y_U = sin(HA_U*pi/180) * cos(Dec_U*pi/180);
     z_U = sin(Dec_U*pi/180);
	 xhor_U = x_U * sin(35*pi/180) - z_U * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_U = y_U; //Horizontal coordinates
     zhor_U = x_U * cos(35*pi/180) + z_U * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_U  = atan2( yhor_U, xhor_U )*180/pi + 180; //local azimuth
     alt_U = atan2( zhor_U, sqrt(xhor_U*xhor_U+yhor_U*yhor_U) )*180/pi; //local altitude
	 alt_topoc_U = alt_U - par_U * cos(alt_U*pi/180); //topocentric altitude
	 gclat_U = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_U   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_U = atan( tan(gclat_U*pi/180) / cos(HA_U*pi/180) ); //an auxiliary angle
	 topRA_U   = RA_U  - par_U * rho_U * cos(gclat_U*pi/180) * sin(HA_U*pi/180) / cos(Dec_U*pi/180);
     topDecl_U = Dec_U - par_U * rho_U * sin(gclat_U*pi/180) * sin(g_U - Dec_U*pi/180) / sin(g_U);
	 
	  UT_U_in_south = ( topRA_U - GMST0 - 51.5 );
	while(UT_U_in_south>360)
		UT_U_in_south = UT_U_in_south - 360;
	while(UT_U_in_south<0)
		UT_U_in_south = UT_U_in_south +360;
	 cosLHA_U = (sin((-0.583)*pi/180) - sin(35.0*pi/180) * sin(topDecl_U*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_U*pi/180));

	if(cosLHA_U<-1.0){
		cin>>D;
	}
	if(cosLHA_U>1.0){
		cin>>D;
	}
	 set_U;
	rise_U;
	if(cosLHA_U<1 && cosLHA_U>-1)
	{
		LHA_U = acos(cosLHA_U)*180/pi; //Time in degrees
		if(LHA_U>(180))
			LHA_U = LHA_U- 360;
		if(LHA_U<-180)
			LHA_U = LHA_U + 360;
		set_U = UT_U_in_south + LHA_U;
	}
	while(set_U>(360))
		set_U = set_U - 360;
	while(set_U<0)
		set_U = set_U + 360;
	set_U = set_U/15.04107  +4.5;
	perturbed_set_time_U = set_U/24;
	}
	perturbed_set_time_U = perturbed_set_time_U*24;

	
	//////////////////////////////////////////////////////////////////////////////////PERTURBING RISE

	for(int rise=0;rise<4;rise++){
	d = reald ; //day
	d = d + perturbed_rise_time_M; //perturbing time
	N_U=74.0005 + 1.3978E-5 * d;  //longitude of the ascending node
    i_U=0.7733 + 1.9E-8 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    w_U=96.6612 + 3.0565E-5 * d;  //argument of perihelion
    a_U=19.18171 - 1.55E-8 * d;  //semi-major axis, or mean distance from Sun
    e_U=0.047318 + 7.45E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M_U=142.5905 + 0.011725806 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_U>(360))
		M_U = M_U - 360;
	while(M_U<0)
		M_U = M_U + 360;
	
	  E_U = M_U + e_U*(180/pi) * sin(M_U*pi/180) * ( 1.0 + e_U * cos(M_U*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_U>0.05){
        double E0_U=E_U;
        double E1_U=E_U;
        while(fabs(E0_U-E1_U)>0.001){
            E0_U=E1_U;
            E1_U = E_U - ( E_U - e_U*(180/pi) * sin(E_U*pi/180) - M_U ) / ( 1 - e_U * cos(E_U*pi/180) );
        }
    }
     xv_U = a_U * (cos(E_U*pi/180) - e_U); //NOT IMPORTANT
     yv_U = a_U * sqrt(1.0 - e_U*e_U) * sin(E_U*pi/180); //NOT IMPORTANT
     v_U = atan2( yv_U, xv_U ); //Mercury's true anomaly
     r_U = sqrt( xv_U*xv_U + yv_U*yv_U ); //Mercury's distance
     xh_U = r_U * ( cos(N_U*pi/180) * cos(v_U+w_U*pi/180) - sin(N_U*pi/180) * sin(v_U+w_U*pi/180) * cos(i_U*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_U = r_U * ( sin(N_U*pi/180) * cos(v_U+w_U*pi/180) + cos(N_U*pi/180) * sin(v_U+w_U*pi/180) * cos(i_U*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_U = r_U * ( sin(v_U+w_U*pi/180) * sin(i_U*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_U = atan2( yh_U, xh_U )*180/pi; //Mercury's ecliptic longitude
     latecl_U = atan2( zh_U, sqrt(xh_U*xh_U+yh_U*yh_U) )*180/pi; //Mercury's ecliptic latitude
	 xg_U = r_U * cos(lonecl_U*pi/180) * cos(latecl_U*pi/180) +xs; //geocentric position
	 yg_U = r_U * sin(lonecl_U*pi/180) * cos(latecl_U*pi/180) +ys; //geocentric position
	 zg_U = r_U * sin(latecl_U*pi/180); //geocentric position
	 xe_U = xg_U; //equatorial coordinates
     ye_U = yg_U * cos(ecl*pi/180) - zg_U * sin(ecl*pi/180); //equatorial coordinates
     ze_U = yg_U * sin(ecl*pi/180) + zg_U * cos(ecl*pi/180); //equatorial coordinates
     RA_U  = atan2( ye_U, xe_U )*180/pi; //Right Ascension
     Dec_U = atan2( ze_U, sqrt(xe_U*xe_U+ye_U*ye_U) )*180/pi; //Declination
	 par_U = (8.794/3600) / r_U;
	 HA_U = LST - RA_U; //hour angle of the Moon
	 x_U = cos(HA_U*pi/180) * cos(Dec_U*pi/180);
     y_U = sin(HA_U*pi/180) * cos(Dec_U*pi/180);
     z_U = sin(Dec_U*pi/180);
	 xhor_U = x_U * sin(35*pi/180) - z_U * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_U = y_U; //Horizontal coordinates
     zhor_U = x_U * cos(35*pi/180) + z_U * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_U  = atan2( yhor_U, xhor_U )*180/pi + 180; //local azimuth
     alt_U = atan2( zhor_U, sqrt(xhor_U*xhor_U+yhor_U*yhor_U) )*180/pi; //local altitude
	 alt_topoc_U = alt_U - par_U * cos(alt_U*pi/180); //topocentric altitude
	 gclat_U = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_U   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_U = atan( tan(gclat_U*pi/180) / cos(HA_U*pi/180) ); //an auxiliary angle
	 topRA_U   = RA_U  - par_U * rho_U * cos(gclat_U*pi/180) * sin(HA_U*pi/180) / cos(Dec_U*pi/180);
     topDecl_U = Dec_U - par_U * rho_U * sin(gclat_U*pi/180) * sin(g_U - Dec_U*pi/180) / sin(g_U);
	
	  UT_U_in_south = ( topRA_U - GMST0 - 51.5 );
	while(UT_U_in_south>360)
		UT_U_in_south = UT_U_in_south - 360;
	while(UT_U_in_south<0)
		UT_U_in_south = UT_U_in_south +360;
	 cosLHA_U = (sin(-0.583*pi/180) - sin(35.0*pi/180) * sin(topDecl_U*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_U*pi/180));

	if(cosLHA_U<-1.0){
		cin>>D;
	}
	if(cosLHA_U>1.0){
		cin>>D;
	}
	if(cosLHA_U<1 && cosLHA_U>-1)
	{
		LHA_U = acos(cosLHA_U)*180/pi; //Time in degrees
		if(LHA_U>(180))
			LHA_U = LHA_U- 360;
		if(LHA_U<-180)
			LHA_U = LHA_U + 360;
		rise_U = UT_U_in_south - LHA_U;
	}
	while(rise_U>(360))
		rise_U = rise_U - 360;
	while(rise_U<0)
		rise_U = rise_U + 360;
	rise_U = rise_U/15.04107  + 4.5;
	perturbed_rise_time_U = rise_U/24;
	}
	perturbed_rise_time_U = perturbed_rise_time_U*24;
	UT_U_in_south = UT_U_in_south/15;
	UT_U_in_south = UT_U_in_south + 4.5;
	double meridianhour_U = (int)UT_U_in_south;
	double meridianminute_U = (int)((UT_U_in_south - meridianhour_U) * 60);
	double meridiansecond_U = (int)((((UT_U_in_south - meridianhour_U) * 60) - meridianminute_U) * 60);
	double sethour_U = (int)set_U;
	double setminute_U = (int)((set_U - sethour_U) * 60);
	double setsecond_U = (int)((((set_U - sethour_U) * 60) - setminute_U) * 60);
	double risehour_U = (int)rise_U;
	double riseminute_U = (int)((rise_U - risehour_U) * 60);
	double risesecond_U = (int)((((rise_U - risehour_U) * 60) - riseminute_U) * 60);
	cout<<"Uranus:"<<endl<< "setting time: "<<sethour_U<<" : "<<setminute_U<<" : "<<setsecond_U<<endl<<"rising time: "<<risehour_U<<" : "<<riseminute_U<<" : "<<risesecond_U<<endl<<"meridian time:"<<meridianhour_U<<" : "<<meridianminute_U<<" : "<<meridiansecond_U<<endl;

	////////////////////////////////Neptune's rise and set times////////////////////////////////
	
		d = reald + 0.3125;
	N_N=131.7806 + 3.0173E-5 * d;  //longitude of the ascending node
    i_N=1.7700 - 2.55E-7 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    w_N=272.8461 - 6.027E-6 * d;  //argument of perihelion
    a_N=30.05826 + 3.313E-8 * d;  //semi-major axis, or mean distance from Sun
    e_N=0.008606 + 2.15E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M_N=260.2471 + 0.005995147 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_N>(360))
		M_N = M_N - 360;
	while(M_N<0)
		M_N = M_N + 360;
    
	  E_N = M_N + e_N*(180/pi) * sin(M_N*pi/180) * ( 1.0 + e_N * cos(M_N*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_N>0.05){
        double E0_N=E_N;
        double E1_N=E_N;
        while(fabs(E0_N-E1_N)>0.001){
            E0_N=E1_N;
            E1_N = E_N - ( E_N - e_N*(180/pi) * sin(E_N*pi/180) - M_N ) / ( 1 - e_N * cos(E_N*pi/180) );
        }
    }
     xv_N = a_N * (cos(E_N*pi/180) - e_N); //NOT IMPORTANT
     yv_N = a_N * sqrt(1.0 - e_N*e_N) * sin(E_N*pi/180); //NOT IMPORTANT
     v_N = atan2( yv_N, xv_N ); //Mercury's true anomaly
     r_N = sqrt( xv_N*xv_N + yv_N*yv_N ); //Mercury's distance
     xh_N = r_N * ( cos(N_N*pi/180) * cos(v_N+w_N*pi/180) - sin(N_N*pi/180) * sin(v_N+w_N*pi/180) * cos(i_N*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_N = r_N * ( sin(N_N*pi/180) * cos(v_N+w_N*pi/180) + cos(N_N*pi/180) * sin(v_N+w_N*pi/180) * cos(i_N*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_N = r_N * ( sin(v_N+w_N*pi/180) * sin(i_N*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_N = atan2( yh_N, xh_N )*180/pi; //Mercury's ecliptic longitude
     latecl_N = atan2( zh_N, sqrt(xh_N*xh_N+yh_N*yh_N) )*180/pi; //Mercury's ecliptic latitude
	 xg_N = r_N * cos(lonecl_N*pi/180) * cos(latecl_N*pi/180) +xs; //geocentric position
	 yg_N = r_N * sin(lonecl_N*pi/180) * cos(latecl_N*pi/180) +ys; //geocentric position
	 zg_N = r_N * sin(latecl_N*pi/180); //geocentric position
	 xe_N = xg_N; //equatorial coordinates
     ye_N = yg_N * cos(ecl*pi/180) - zg_N * sin(ecl*pi/180); //equatorial coordinates
     ze_N = yg_N * sin(ecl*pi/180) + zg_N * cos(ecl*pi/180); //equatorial coordinates
     RA_N  = atan2( ye_N, xe_N )*180/pi; //Right Ascension
     Dec_N = atan2( ze_N, sqrt(xe_N*xe_N+ye_N*ye_N) )*180/pi; //Declination
	 par_N = (8.794/3600) / r_N;
	 HA_N = LST - RA_N; //hour angle of the Moon
	 x_N = cos(HA_N*pi/180) * cos(Dec_N*pi/180);
     y_N = sin(HA_N*pi/180) * cos(Dec_N*pi/180);
     z_N = sin(Dec_N*pi/180);
	 xhor_N = x_N * sin(35*pi/180) - z_N * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_N = y_N; //Horizontal coordinates
     zhor_N = x_N * cos(35*pi/180) + z_N * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_N  = atan2( yhor_N, xhor_N )*180/pi + 180; //local azimuth
     alt_N = atan2( zhor_N, sqrt(xhor_N*xhor_N+yhor_N*yhor_N) )*180/pi; //local altitude
	 alt_topoc_N = alt_N - par_N * cos(alt_N*pi/180); //topocentric altitude
	 gclat_N = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_N   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_N = atan( tan(gclat_N*pi/180) / cos(HA_N*pi/180) ); //an auxiliary angle
	 topRA_N   = RA_N  - par_N * rho_N * cos(gclat_N*pi/180) * sin(HA_N*pi/180) / cos(Dec_N*pi/180);
     topDecl_N = Dec_N - par_N * rho_N * sin(gclat_N*pi/180) * sin(g_N - Dec_N*pi/180) / sin(g_N);
	
	 double UT_N_in_south = ( topRA_N - GMST0 - 51.5 );
	while(UT_N_in_south>360)
		UT_N_in_south = UT_N_in_south - 360;
	while(UT_N_in_south<0)
		UT_N_in_south = UT_N_in_south +360;
	double cosLHA_N = (sin((-0.583)*pi/180) - sin(35.0*pi/180) * sin(topDecl_N*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_N*pi/180));
	double LHA_N;
	if(cosLHA_N<-1.0){
		cin>>D;
	}
	if(cosLHA_N>1.0){
		cin>>D;
	}
	double set_N;
	double rise_N;
	if(cosLHA_N<1 && cosLHA_N>-1)
	{
		LHA_N = acos(cosLHA_N)*180/pi; //Time in degrees
		if(LHA_N>(180))
			LHA_N = LHA_N- 360;
		if(LHA_N<-180)
			LHA_N = LHA_N + 360;
		set_N = UT_N_in_south + LHA_N;
		rise_N = UT_N_in_south - LHA_N;
	}
	while(set_N>(360))
		set_N = set_N - 360;
	while(set_N<0)
		set_N = set_N + 360;
	while(rise_N>(360))
		rise_N = rise_N - 360;
	while(rise_N<0)
		rise_N = rise_N + 360;
	set_N = set_N/15.04107  +4.5;
	rise_N = rise_N/15.04107  +4.5;
	//////////////////////////////////////////////////////////////////////////////////PERTURBING SET
	double perturbed_rise_time_N = rise_N/24;
	double perturbed_set_time_N = set_N/24;

	for(int set=0;set<4;set++){
	d = reald ; //day
	d = d + perturbed_set_time_N; //perturbing time
	N_N=131.7806 + 3.0173E-5 * d;  //longitude of the ascending node
    i_N=1.7700 - 2.55E-7 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    w_N=272.8461 - 6.027E-6 * d;  //argument of perihelion
    a_N=30.05826 + 3.313E-8 * d;  //semi-major axis, or mean distance from Sun
    e_N=0.008606 + 2.15E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M_N=260.2471 + 0.005995147 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_N>(360))
		M_N = M_N - 360;
	while(M_N<0)
		M_N = M_N + 360;
    
	  E_N = M_N + e_N*(180/pi) * sin(M_N*pi/180) * ( 1.0 + e_N * cos(M_N*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_N>0.05){
        double E0_N=E_N;
        double E1_N=E_N;
        while(fabs(E0_N-E1_N)>0.001){
            E0_N=E1_N;
            E1_N = E_N - ( E_N - e_N*(180/pi) * sin(E_N*pi/180) - M_N ) / ( 1 - e_N * cos(E_N*pi/180) );
        }
    }
     xv_N = a_N * (cos(E_N*pi/180) - e_N); //NOT IMPORTANT
     yv_N = a_N * sqrt(1.0 - e_N*e_N) * sin(E_N*pi/180); //NOT IMPORTANT
     v_N = atan2( yv_N, xv_N ); //Mercury's true anomaly
     r_N = sqrt( xv_N*xv_N + yv_N*yv_N ); //Mercury's distance
     xh_N = r_N * ( cos(N_N*pi/180) * cos(v_N+w_N*pi/180) - sin(N_N*pi/180) * sin(v_N+w_N*pi/180) * cos(i_N*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_N = r_N * ( sin(N_N*pi/180) * cos(v_N+w_N*pi/180) + cos(N_N*pi/180) * sin(v_N+w_N*pi/180) * cos(i_N*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_N = r_N * ( sin(v_N+w_N*pi/180) * sin(i_N*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_N = atan2( yh_N, xh_N )*180/pi; //Mercury's ecliptic longitude
     latecl_N = atan2( zh_N, sqrt(xh_N*xh_N+yh_N*yh_N) )*180/pi; //Mercury's ecliptic latitude
	 xg_N = r_N * cos(lonecl_N*pi/180) * cos(latecl_N*pi/180) +xs; //geocentric position
	 yg_N = r_N * sin(lonecl_N*pi/180) * cos(latecl_N*pi/180) +ys; //geocentric position
	 zg_N = r_N * sin(latecl_N*pi/180); //geocentric position
	 xe_N = xg_N; //equatorial coordinates
     ye_N = yg_N * cos(ecl*pi/180) - zg_N * sin(ecl*pi/180); //equatorial coordinates
     ze_N = yg_N * sin(ecl*pi/180) + zg_N * cos(ecl*pi/180); //equatorial coordinates
     RA_N  = atan2( ye_N, xe_N )*180/pi; //Right Ascension
     Dec_N = atan2( ze_N, sqrt(xe_N*xe_N+ye_N*ye_N) )*180/pi; //Declination
	 par_N = (8.794/3600) / r_N;
	 HA_N = LST - RA_N; //hour angle of the Moon
	 x_N = cos(HA_N*pi/180) * cos(Dec_N*pi/180);
     y_N = sin(HA_N*pi/180) * cos(Dec_N*pi/180);
     z_N = sin(Dec_N*pi/180);
	 xhor_N = x_N * sin(35*pi/180) - z_N * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_N = y_N; //Horizontal coordinates
     zhor_N = x_N * cos(35*pi/180) + z_N * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_N  = atan2( yhor_N, xhor_N )*180/pi + 180; //local azimuth
     alt_N = atan2( zhor_N, sqrt(xhor_N*xhor_N+yhor_N*yhor_N) )*180/pi; //local altitude
	 alt_topoc_N = alt_N - par_N * cos(alt_N*pi/180); //topocentric altitude
	 gclat_N = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_N   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_N = atan( tan(gclat_N*pi/180) / cos(HA_N*pi/180) ); //an auxiliary angle
	 topRA_N   = RA_N  - par_N * rho_N * cos(gclat_N*pi/180) * sin(HA_N*pi/180) / cos(Dec_N*pi/180);
     topDecl_N = Dec_N - par_N * rho_N * sin(gclat_N*pi/180) * sin(g_N - Dec_N*pi/180) / sin(g_N);
	 
	  UT_N_in_south = ( topRA_N - GMST0 - 51.5 );
	while(UT_N_in_south>360)
		UT_N_in_south = UT_N_in_south - 360;
	while(UT_N_in_south<0)
		UT_N_in_south = UT_N_in_south +360;
	 cosLHA_N = (sin((-0.583)*pi/180) - sin(35.0*pi/180) * sin(topDecl_N*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_N*pi/180));

	if(cosLHA_N<-1.0){
		cin>>D;
	}
	if(cosLHA_N>1.0){
		cin>>D;
	}
	 set_N;
	rise_N;
	if(cosLHA_N<1 && cosLHA_N>-1)
	{
		LHA_N = acos(cosLHA_N)*180/pi; //Time in degrees
		if(LHA_N>(180))
			LHA_N = LHA_N- 360;
		if(LHA_N<-180)
			LHA_N = LHA_N + 360;
		set_N = UT_N_in_south + LHA_N;
	}
	while(set_N>(360))
		set_N = set_N - 360;
	while(set_N<0)
		set_N = set_N + 360;
	set_N = set_N/15.04107  +4.5;
	perturbed_set_time_N = set_N/24;
	}
	perturbed_set_time_N = perturbed_set_time_N*24;

	
	//////////////////////////////////////////////////////////////////////////////////PERTURBING RISE

	for(int rise=0;rise<4;rise++){
	d = reald ; //day
	d = d + perturbed_rise_time_M; //perturbing time
	N_N=131.7806 + 3.0173E-5 * d;  //longitude of the ascending node
    i_N=1.7700 - 2.55E-7 * d;  //inclination to the ecliptic (plane of the Earth's orbit)
    w_N=272.8461 - 6.027E-6 * d;  //argument of perihelion
    a_N=30.05826 + 3.313E-8 * d;  //semi-major axis, or mean distance from Sun
    e_N=0.008606 + 2.15E-9 * d;  //eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M_N=260.2471 + 0.005995147 * d;  //mean anomaly (0 at perihelion; increases uniformly with time)
		while(M_N>(360))
		M_N = M_N - 360;
	while(M_N<0)
		M_N = M_N + 360;
    
	  E_N = M_N + e_N*(180/pi) * sin(M_N*pi/180) * ( 1.0 + e_N * cos(M_N*pi/180) );  //Mercury's eccentric anomaly E_S
    if(e_N>0.05){
        double E0_N=E_N;
        double E1_N=E_N;
        while(fabs(E0_N-E1_N)>0.001){
            E0_N=E1_N;
            E1_N = E_N - ( E_N - e_N*(180/pi) * sin(E_N*pi/180) - M_N ) / ( 1 - e_N * cos(E_N*pi/180) );
        }
    }
     xv_N = a_N * (cos(E_N*pi/180) - e_N); //NOT IMPORTANT
     yv_N = a_N * sqrt(1.0 - e_N*e_N) * sin(E_N*pi/180); //NOT IMPORTANT
     v_N = atan2( yv_N, xv_N ); //Mercury's true anomaly
     r_N = sqrt( xv_N*xv_N + yv_N*yv_N ); //Mercury's distance
     xh_N = r_N * ( cos(N_N*pi/180) * cos(v_N+w_N*pi/180) - sin(N_N*pi/180) * sin(v_N+w_N*pi/180) * cos(i_N*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     yh_N = r_N * ( sin(N_N*pi/180) * cos(v_N+w_N*pi/180) + cos(N_N*pi/180) * sin(v_N+w_N*pi/180) * cos(i_N*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
     zh_N = r_N * ( sin(v_N+w_N*pi/180) * sin(i_N*pi/180) ); //Mercury's heliocentric (Sun-centered) position in the ecliptic coordinate system
	 lonecl_N = atan2( yh_N, xh_N )*180/pi; //Mercury's ecliptic longitude
     latecl_N = atan2( zh_N, sqrt(xh_N*xh_N+yh_N*yh_N) )*180/pi; //Mercury's ecliptic latitude
	 xg_N = r_N * cos(lonecl_N*pi/180) * cos(latecl_N*pi/180) +xs; //geocentric position
	 yg_N = r_N * sin(lonecl_N*pi/180) * cos(latecl_N*pi/180) +ys; //geocentric position
	 zg_N = r_N * sin(latecl_N*pi/180); //geocentric position
	 xe_N = xg_N; //equatorial coordinates
     ye_N = yg_N * cos(ecl*pi/180) - zg_N * sin(ecl*pi/180); //equatorial coordinates
     ze_N = yg_N * sin(ecl*pi/180) + zg_N * cos(ecl*pi/180); //equatorial coordinates
     RA_N  = atan2( ye_N, xe_N )*180/pi; //Right Ascension
     Dec_N = atan2( ze_N, sqrt(xe_N*xe_N+ye_N*ye_N) )*180/pi; //Declination
	 par_N = (8.794/3600) / r_N;
	 HA_N = LST - RA_N; //hour angle of the Moon
	 x_N = cos(HA_N*pi/180) * cos(Dec_N*pi/180);
     y_N = sin(HA_N*pi/180) * cos(Dec_N*pi/180);
     z_N = sin(Dec_N*pi/180);
	 xhor_N = x_N * sin(35*pi/180) - z_N * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_N = y_N; //Horizontal coordinates
     zhor_N = x_N * cos(35*pi/180) + z_N * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_N  = atan2( yhor_N, xhor_N )*180/pi + 180; //local azimuth
     alt_N = atan2( zhor_N, sqrt(xhor_N*xhor_N+yhor_N*yhor_N) )*180/pi; //local altitude
	 alt_topoc_N = alt_N - par_N * cos(alt_N*pi/180); //topocentric altitude
	 gclat_N = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_N   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_N = atan( tan(gclat_N*pi/180) / cos(HA_N*pi/180) ); //an auxiliary angle
	 topRA_N   = RA_N  - par_N * rho_N * cos(gclat_N*pi/180) * sin(HA_N*pi/180) / cos(Dec_N*pi/180);
     topDecl_N = Dec_N - par_N * rho_N * sin(gclat_N*pi/180) * sin(g_N - Dec_N*pi/180) / sin(g_N);
	
	  UT_N_in_south = ( topRA_N - GMST0 - 51.5 );
	while(UT_N_in_south>360)
		UT_N_in_south = UT_N_in_south - 360;
	while(UT_N_in_south<0)
		UT_N_in_south = UT_N_in_south +360;
	 cosLHA_N = (sin(-0.583*pi/180) - sin(35.0*pi/180) * sin(topDecl_N*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_N*pi/180));

	if(cosLHA_N<-1.0){
		cin>>D;
	}
	if(cosLHA_N>1.0){
		cin>>D;
	}
	if(cosLHA_N<1 && cosLHA_N>-1)
	{
		LHA_N = acos(cosLHA_N)*180/pi; //Time in degrees
		if(LHA_N>(180))
			LHA_N = LHA_N- 360;
		if(LHA_N<-180)
			LHA_N = LHA_N + 360;
		rise_N = UT_N_in_south - LHA_N;
	}
	while(rise_N>(360))
		rise_N = rise_N - 360;
	while(rise_N<0)
		rise_N = rise_N + 360;
	rise_N = rise_N/15.04107  + 4.5;
	perturbed_rise_time_N = rise_N/24;
	}
	perturbed_rise_time_N = perturbed_rise_time_N*24;
	UT_N_in_south = UT_N_in_south/15;
	UT_N_in_south = UT_N_in_south + 4.5;
	double meridianhour_N = (int)UT_N_in_south;
	double meridianminute_N = (int)((UT_N_in_south - meridianhour_N) * 60);
	double meridiansecond_N = (int)((((UT_N_in_south - meridianhour_N) * 60) - meridianminute_N) * 60);
	double sethour_N = (int)set_N;
	double setminute_N = (int)((set_N - sethour_N) * 60);
	double setsecond_N = (int)((((set_N - sethour_N) * 60) - setminute_N) * 60);
	double risehour_N = (int)rise_N;
	double riseminute_N = (int)((rise_N - risehour_N) * 60);
	double risesecond_N = (int)((((rise_N - risehour_N) * 60) - riseminute_N) * 60);
	cout<<"Neptune:"<<endl<< "setting time: "<<sethour_N<<" : "<<setminute_N<<" : "<<setsecond_N<<endl<<"rising time: "<<risehour_N<<" : "<<riseminute_N<<" : "<<risesecond_N<<endl<<"meridian time:"<<meridianhour_N<<" : "<<meridianminute_N<<" : "<<meridiansecond_N<<endl;

	////////////////////////////////Pluto's rise and set times////////////////////////////////
	d = reald + 0.3125;
	S  =   50.03  +  0.033459652 * d;
     P  =  238.95  +  0.003968789 * d;
	 lonecl_P = 238.9508  +  0.00400703 * d
            - 19.799 * sin(P*pi/180)     + 19.848 * cos(P*pi/180)
             + 0.897 * sin(2*P*pi/180)    - 4.956 * cos(2*P*pi/180)
             + 0.610 * sin(3*P*pi/180)    + 1.211 * cos(3*P*pi/180)
             - 0.341 * sin(4*P*pi/180)    - 0.190 * cos(4*P*pi/180)
             + 0.128 * sin(5*P*pi/180)    - 0.034 * cos(5*P*pi/180)
             - 0.038 * sin(6*P*pi/180)    + 0.031 * cos(6*P*pi/180)
             + 0.020 * sin(S*pi/180-P*pi/180)    - 0.010 * cos(S*pi/180-P*pi/180);
	  latecl_P = 453 * sin(P*pi/180)     - 14.975 * cos(P*pi/180)
             + 3.527 * sin(2*P*pi/180)    + 1.673 * cos(2*P*pi/180)
             - 1.051 * sin(3*P*pi/180)    + 0.328 * cos(3*P*pi/180)
             + 0.179 * sin(4*P*pi/180)    - 0.292 * cos(4*P*pi/180)
             + 0.019 * sin(5*P*pi/180)    + 0.100 * cos(5*P*pi/180)
             - 0.031 * sin(6*P*pi/180)    - 0.026 * cos(6*P*pi/180)
                                   + 0.011 * cos(S*pi/180-P*pi/180);
	  r_P     =  40.72
           + 6.68 * sin(P*pi/180)       + 6.90 * cos(P*pi/180)
           - 1.18 * sin(2*P*pi/180)     - 0.03 * cos(2*P*pi/180)
           + 0.15 * sin(3*P*pi/180)     - 0.14 * cos(3*P*pi/180);
	 
	while(lonecl_P>90)
		lonecl_P = lonecl_P - 90;
	while(lonecl_P<-90)
		lonecl_P = lonecl_P + 90;

	while(latecl_P>90)
		latecl_P = latecl_P - 90;
	while(latecl_P<-90)
		latecl_P = latecl_P + 90;
     xg_P = r_P * cos(lonecl_P*pi/180) * cos(latecl_P*pi/180) +xs; //geocentric position
	 yg_P = r_P * sin(lonecl_P*pi/180) * cos(latecl_P*pi/180) +ys; //geocentric position
	 zg_P = r_P * sin(latecl_P*pi/180); //geocentric position
	 xe_P = xg_P; //equatorial coordinates
     ye_P = yg_P * cos(ecl*pi/180) - zg_P * sin(ecl*pi/180); //equatorial coordinates
     ze_P = yg_P * sin(ecl*pi/180) + zg_P * cos(ecl*pi/180); //equatorial coordinates
     RA_P  = atan2( ye_P, xe_P )*180/pi; //Right Ascension
     Dec_P = atan2( ze_P, sqrt(xe_P*xe_P+ye_P*ye_P) )*180/pi; //Declination
	 par_P = (8.794/3600) / r_P;
	 HA_P = LST - RA_P; //hour angle of the Moon
	 x_P = cos(HA_P*pi/180) * cos(Dec_P*pi/180);
     y_P = sin(HA_P*pi/180) * cos(Dec_P*pi/180);
     z_P = sin(Dec_P*pi/180);
	 xhor_P = x_P * sin(35*pi/180) - z_P * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_P = y_P; //Horizontal coordinates
     zhor_P = x_P * cos(35*pi/180) + z_P * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_P  = atan2( yhor_P, xhor_P )*180/pi + 180; //local azimuth
     alt_P = atan2( zhor_P, sqrt(xhor_P*xhor_P+yhor_P*yhor_P) )*180/pi; //local altitude
	 alt_topoc_P = alt_P - par_P * cos(alt_P*pi/180); //topocentric altitude
	 gclat_P = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_P   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_P = atan( tan(gclat_P*pi/180) / cos(HA_P*pi/180) ); //an auxiliary angle
	 topRA_P   = RA_P  - par_P * rho_P * cos(gclat_P*pi/180) * sin(HA_P*pi/180) / cos(Dec_P*pi/180);
     topDecl_P = Dec_P - par_P * rho_P * sin(gclat_P*pi/180) * sin(g_P - Dec_P*pi/180) / sin(g_P);
	
	 double UT_P_in_south = ( topRA_P - GMST0 - 51.5 );
	while(UT_P_in_south>360)
		UT_P_in_south = UT_P_in_south - 360;
	while(UT_P_in_south<0)
		UT_P_in_south = UT_P_in_south +360;
	double cosLHA_P = (sin((-0.583)*pi/180) - sin(35.0*pi/180) * sin(topDecl_P*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_P*pi/180));
	double LHA_P;
	if(cosLHA_P<-1.0){
		cin>>D;
	}
	if(cosLHA_P>1.0){
		cin>>D;
	}
	double set_P;
	double rise_P;
	if(cosLHA_P<1 && cosLHA_P>-1)
	{
		LHA_P = acos(cosLHA_P)*180/pi; //Time in degrees
		if(LHA_P>(180))
			LHA_P = LHA_P- 360;
		if(LHA_P<-180)
			LHA_P = LHA_P + 360;
		set_P = UT_P_in_south + LHA_P;
		rise_P = UT_P_in_south - LHA_P;
	}
	while(set_P>(360))
		set_P = set_P - 360;
	while(set_P<0)
		set_P = set_P + 360;
	while(rise_P>(360))
		rise_P = rise_P - 360;
	while(rise_P<0)
		rise_P = rise_P + 360;
	set_P = set_P/15.04107  +4.5;
	rise_P = rise_P/15.04107  +4.5;
	//////////////////////////////////////////////////////////////////////////////////PERTURBING SET
	double perturbed_rise_time_P = rise_P/24;
	double perturbed_set_time_P = set_P/24;

	for(int set=0;set<4;set++){
	d = reald ; //day
	d = d + perturbed_set_time_P; //perturbing time
	S  =   50.03  +  0.033459652 * d;
     P  =  238.95  +  0.003968789 * d;
	 lonecl_P = 238.9508  +  0.00400703 * d
            - 19.799 * sin(P*pi/180)     + 19.848 * cos(P*pi/180)
             + 0.897 * sin(2*P*pi/180)    - 4.956 * cos(2*P*pi/180)
             + 0.610 * sin(3*P*pi/180)    + 1.211 * cos(3*P*pi/180)
             - 0.341 * sin(4*P*pi/180)    - 0.190 * cos(4*P*pi/180)
             + 0.128 * sin(5*P*pi/180)    - 0.034 * cos(5*P*pi/180)
             - 0.038 * sin(6*P*pi/180)    + 0.031 * cos(6*P*pi/180)
             + 0.020 * sin(S*pi/180-P*pi/180)    - 0.010 * cos(S*pi/180-P*pi/180);
	  latecl_P = 453 * sin(P*pi/180)     - 14.975 * cos(P*pi/180)
             + 3.527 * sin(2*P*pi/180)    + 1.673 * cos(2*P*pi/180)
             - 1.051 * sin(3*P*pi/180)    + 0.328 * cos(3*P*pi/180)
             + 0.179 * sin(4*P*pi/180)    - 0.292 * cos(4*P*pi/180)
             + 0.019 * sin(5*P*pi/180)    + 0.100 * cos(5*P*pi/180)
             - 0.031 * sin(6*P*pi/180)    - 0.026 * cos(6*P*pi/180)
                                   + 0.011 * cos(S*pi/180-P*pi/180);
	  r_P     =  40.72
           + 6.68 * sin(P*pi/180)       + 6.90 * cos(P*pi/180)
           - 1.18 * sin(2*P*pi/180)     - 0.03 * cos(2*P*pi/180)
           + 0.15 * sin(3*P*pi/180)     - 0.14 * cos(3*P*pi/180);
	 
	while(lonecl_P>90)
		lonecl_P = lonecl_P - 90;
	while(lonecl_P<-90)
		lonecl_P = lonecl_P + 90;

	while(latecl_P>90)
		latecl_P = latecl_P - 90;
	while(latecl_P<-90)
		latecl_P = latecl_P + 90;
     xg_P = r_P * cos(lonecl_P*pi/180) * cos(latecl_P*pi/180) +xs; //geocentric position
	 yg_P = r_P * sin(lonecl_P*pi/180) * cos(latecl_P*pi/180) +ys; //geocentric position
	 zg_P = r_P * sin(latecl_P*pi/180); //geocentric position
	 xe_P = xg_P; //equatorial coordinates
     ye_P = yg_P * cos(ecl*pi/180) - zg_P * sin(ecl*pi/180); //equatorial coordinates
     ze_P = yg_P * sin(ecl*pi/180) + zg_P * cos(ecl*pi/180); //equatorial coordinates
     RA_P  = atan2( ye_P, xe_P )*180/pi; //Right Ascension
     Dec_P = atan2( ze_P, sqrt(xe_P*xe_P+ye_P*ye_P) )*180/pi; //Declination
	 par_P = (8.794/3600) / r_P;
	 HA_P = LST - RA_P; //hour angle of the Moon
	 x_P = cos(HA_P*pi/180) * cos(Dec_P*pi/180);
     y_P = sin(HA_P*pi/180) * cos(Dec_P*pi/180);
     z_P = sin(Dec_P*pi/180);
	 xhor_P = x_P * sin(35*pi/180) - z_P * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_P = y_P; //Horizontal coordinates
     zhor_P = x_P * cos(35*pi/180) + z_P * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_P  = atan2( yhor_P, xhor_P )*180/pi + 180; //local azimuth
     alt_P = atan2( zhor_P, sqrt(xhor_P*xhor_P+yhor_P*yhor_P) )*180/pi; //local altitude
	 (8.794/3600) / r_P;
	 alt_topoc_P = alt_P - par_P * cos(alt_P*pi/180); //topocentric altitude
	 gclat_P = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_P   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_P = atan( tan(gclat_P*pi/180) / cos(HA_P*pi/180) ); //an auxiliary angle
	 topRA_P   = RA_P  - par_P * rho_P * cos(gclat_P*pi/180) * sin(HA_P*pi/180) / cos(Dec_P*pi/180);
     topDecl_P = Dec_P - par_P * rho_P * sin(gclat_P*pi/180) * sin(g_P - Dec_P*pi/180) / sin(g_P);
	

 
	  UT_P_in_south = ( topRA_P - GMST0 - 51.5 );
	while(UT_P_in_south>360)
		UT_P_in_south = UT_P_in_south - 360;
	while(UT_P_in_south<0)
		UT_P_in_south = UT_P_in_south +360;
	 cosLHA_P = (sin((-0.583)*pi/180) - sin(35.0*pi/180) * sin(topDecl_P*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_P*pi/180));

	if(cosLHA_P<-1.0){
		cin>>D;
	}
	if(cosLHA_P>1.0){
		cin>>D;
	}
	 set_P;
	rise_P;
	if(cosLHA_P<1 && cosLHA_P>-1)
	{
		LHA_P = acos(cosLHA_P)*180/pi; //Time in degrees
		if(LHA_P>(180))
			LHA_P = LHA_P- 360;
		if(LHA_P<-180)
			LHA_P = LHA_P + 360;
		set_P = UT_P_in_south + LHA_P;
	}
	while(set_P>(360))
		set_P = set_P - 360;
	while(set_P<0)
		set_P = set_P + 360;
	set_P = set_P/15.04107  +4.5;
	perturbed_set_time_P = set_P/24;
	}
	perturbed_set_time_P = perturbed_set_time_P*24;

	
	//////////////////////////////////////////////////////////////////////////////////PERTURBING RISE

	for(int rise=0;rise<4;rise++){
	d = reald ; //day
	d = d + perturbed_rise_time_M; //perturbing time
    S  =   50.03  +  0.033459652 * d;
     P  =  238.95  +  0.003968789 * d;
	 lonecl_P = 238.9508  +  0.00400703 * d
            - 19.799 * sin(P*pi/180)     + 19.848 * cos(P*pi/180)
             + 0.897 * sin(2*P*pi/180)    - 4.956 * cos(2*P*pi/180)
             + 0.610 * sin(3*P*pi/180)    + 1.211 * cos(3*P*pi/180)
             - 0.341 * sin(4*P*pi/180)    - 0.190 * cos(4*P*pi/180)
             + 0.128 * sin(5*P*pi/180)    - 0.034 * cos(5*P*pi/180)
             - 0.038 * sin(6*P*pi/180)    + 0.031 * cos(6*P*pi/180)
             + 0.020 * sin(S*pi/180-P*pi/180)    - 0.010 * cos(S*pi/180-P*pi/180);
	  latecl_P = 453 * sin(P*pi/180)     - 14.975 * cos(P*pi/180)
             + 3.527 * sin(2*P*pi/180)    + 1.673 * cos(2*P*pi/180)
             - 1.051 * sin(3*P*pi/180)    + 0.328 * cos(3*P*pi/180)
             + 0.179 * sin(4*P*pi/180)    - 0.292 * cos(4*P*pi/180)
             + 0.019 * sin(5*P*pi/180)    + 0.100 * cos(5*P*pi/180)
             - 0.031 * sin(6*P*pi/180)    - 0.026 * cos(6*P*pi/180)
                                   + 0.011 * cos(S*pi/180-P*pi/180);
	  r_P     =  40.72
           + 6.68 * sin(P*pi/180)       + 6.90 * cos(P*pi/180)
           - 1.18 * sin(2*P*pi/180)     - 0.03 * cos(2*P*pi/180)
           + 0.15 * sin(3*P*pi/180)     - 0.14 * cos(3*P*pi/180);
	 
	while(lonecl_P>90)
		lonecl_P = lonecl_P - 90;
	while(lonecl_P<-90)
		lonecl_P = lonecl_P + 90;

	while(latecl_P>90)
		latecl_P = latecl_P - 90;
	while(latecl_P<-90)
		latecl_P = latecl_P + 90;
     xg_P = r_P * cos(lonecl_P*pi/180) * cos(latecl_P*pi/180) +xs; //geocentric position
	 yg_P = r_P * sin(lonecl_P*pi/180) * cos(latecl_P*pi/180) +ys; //geocentric position
	 zg_P = r_P * sin(latecl_P*pi/180); //geocentric position
	 xe_P = xg_P; //equatorial coordinates
     ye_P = yg_P * cos(ecl*pi/180) - zg_P * sin(ecl*pi/180); //equatorial coordinates
     ze_P = yg_P * sin(ecl*pi/180) + zg_P * cos(ecl*pi/180); //equatorial coordinates
     RA_P  = atan2( ye_P, xe_P )*180/pi; //Right Ascension
     Dec_P = atan2( ze_P, sqrt(xe_P*xe_P+ye_P*ye_P) )*180/pi; //Declination
	 par_P = (8.794/3600) / r_P;
	 HA_P = LST - RA_P; //hour angle of the Moon
	 x_P = cos(HA_P*pi/180) * cos(Dec_P*pi/180);
     y_P = sin(HA_P*pi/180) * cos(Dec_P*pi/180);
     z_P = sin(Dec_P*pi/180);
	 xhor_P = x_P * sin(35*pi/180) - z_P * cos(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     yhor_P = y_P; //Horizontal coordinates
     zhor_P = x_P * cos(35*pi/180) + z_P * sin(35*pi/180); //Horizontal coordinates. replace 35 by your local latitude
     az_P  = atan2( yhor_P, xhor_P )*180/pi + 180; //local azimuth
     alt_P = atan2( zhor_P, sqrt(xhor_P*xhor_P+yhor_P*yhor_P) )*180/pi; //local altitude
	 (8.794/3600) / r_P;
	 alt_topoc_P = alt_P - par_P * cos(alt_P*pi/180); //topocentric altitude
	 gclat_P = 35 - 0.1924 * sin(2*35*pi/180); //geocentric latitude (replace 35 by your local latitude)
     rho_P   = 0.99833 + 0.00167 * cos(2*35*pi/180); //distance from the center of the Earth (replace 35 by your local latitude)
	 g_P = atan( tan(gclat_P*pi/180) / cos(HA_P*pi/180) ); //an auxiliary angle
	 topRA_P   = RA_P  - par_P * rho_P * cos(gclat_P*pi/180) * sin(HA_P*pi/180) / cos(Dec_P*pi/180);
     topDecl_P = Dec_P - par_P * rho_P * sin(gclat_P*pi/180) * sin(g_P - Dec_P*pi/180) / sin(g_P);
	

	  UT_P_in_south = ( topRA_P - GMST0 - 51.5 );
	while(UT_P_in_south>360)
		UT_P_in_south = UT_P_in_south - 360;
	while(UT_P_in_south<0)
		UT_P_in_south = UT_P_in_south +360;
	 cosLHA_P = (sin(-0.583*pi/180) - sin(35.0*pi/180) * sin(topDecl_P*pi/180)) / (cos(35.0*pi/180) * cos(topDecl_P*pi/180));

	if(cosLHA_P<-1.0){
		cin>>D;
	}
	if(cosLHA_P>1.0){
		cin>>D;
	}
	if(cosLHA_P<1 && cosLHA_P>-1)
	{
		LHA_P = acos(cosLHA_P)*180/pi; //Time in degrees
		if(LHA_P>(180))
			LHA_P = LHA_P- 360;
		if(LHA_P<-180)
			LHA_P = LHA_P + 360;
		rise_P = UT_P_in_south - LHA_P;
	}
	while(rise_P>(360))
		rise_P = rise_P - 360;
	while(rise_P<0)
		rise_P = rise_P + 360;
	rise_P = rise_P/15.04107  + 4.5;
	perturbed_rise_time_P = rise_P/24;
	}
	perturbed_rise_time_P = perturbed_rise_time_P*24;
	UT_P_in_south = UT_P_in_south/15;
	UT_P_in_south = UT_P_in_south + 4.5;
	double meridianhour_P = (int)UT_P_in_south;
	double meridianminute_P = (int)((UT_P_in_south - meridianhour_P) * 60);
	double meridiansecond_P = (int)((((UT_P_in_south - meridianhour_P) * 60) - meridianminute_P) * 60);
	double sethour_P = (int)set_P;
	double setminute_P = (int)((set_P - sethour_P) * 60);
	double setsecond_P = (int)((((set_P - sethour_P) * 60) - setminute_P) * 60);
	double risehour_P = (int)rise_P;
	double riseminute_P = (int)((rise_P - risehour_P) * 60);
	double risesecond_P = (int)((((rise_P - risehour_P) * 60) - riseminute_P) * 60);
	cout<<"Pluto:"<<endl<< "setting time: "<<sethour_P<<" : "<<setminute_P<<" : "<<setsecond_P<<endl<<"rising time: "<<risehour_P<<" : "<<riseminute_P<<" : "<<risesecond_P<<endl<<"meridian time:"<<meridianhour_P<<" : "<<meridianminute_P<<" : "<<meridiansecond_P<<endl;
	
    ////////////////////////////////////////////////////////////////



	cin>>D;
    return 0;
}