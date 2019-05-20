#include <iostream>
#include <ginac/ginac.h>
#include "wrpfinite.h"
#include <fstream>
#include <string>

using namespace std;
using namespace GiNaC;

extern "C" double wrpfinite(double y, double z, double logm, int npart){

	string textread,fline;

	ex xx=1.0-y-z;
	ex yy=y;
	ex zz=z;
	ex logmhmu=logm;
	symtab table;
	table["x"]=xx;
	table["y"]=yy;
	table["z"]=zz;
	table["logmh2omu2"]=logmhmu;

    table["Gzy"]=G(lst(-zz),yy);
    table["Gomzy"]=-G(lst(1 - zz),yy);
    table["Gomxy"]=-G(lst(1 - xx),yy);
    table["Gomyx"]=-G(lst(1 - yy),xx);
    table["Gomzx"]=-G(lst(1 - zz),xx);
    table["Gomyz"]=-G(lst(1 - yy),zz);
    table["Gomxz"]=-G(lst(1 - xx),zz);

    table["G0x"]=G(lst(0),xx);
    table["G1x"]=-G(lst(1),xx);
    table["G00x"]=G(lst(0,0),xx);
    table["G10x"]=-G(lst(1,0),xx);
    table["G01x"]=-G(lst(0,1),xx);
    table["G11x"]=G(lst(1,1),xx);
    table["G0omyx"]=-G(lst(0,1 - yy),xx);
    table["Gomy0x"]=-G(lst(1 - yy,0),xx);
    table["Gomyomyx"]=G(lst(1 - yy,1 - yy),xx);
    table["Gyomyx"]=-G(lst(-yy,1 - yy),xx);
    table["G0omzx"]=-G(lst(0,1 - zz),xx);
    table["Gomz0x"]=-G(lst(1 - zz,0),xx);
    table["Gomzomzx"]=G(lst(1 - zz,1 - zz),xx);
    table["Gzomzx"]=-G(lst(-zz,1 - zz),xx);
    table["G000x"]=G(lst(0,0,0),xx);
    table["G010x"]=-G(lst(0,1,0),xx);
    table["G100x"]=-G(lst(1,0,0),xx);
    table["G110x"]=G(lst(1,1,0),xx);
    table["Gomy10x"]=G(lst(1 - yy,1,0),xx);
    table["Gomz10x"]=G(lst(1 - zz,1,0),xx);
    table["G0010x"]=-G(lst(0,0,1,0),xx);
    table["G0100x"]=-G(lst(0,1,0,0),xx);
    table["G0110x"]=G(lst(0,1,1,0),xx);
    table["G1000x"]=-G(lst(1,0,0,0),xx);
    table["G1010x"]=G(lst(1,0,1,0),xx);
    table["G1100x"]=G(lst(1,1,0,0),xx);
    table["G1110x"]=-G(lst(1,1,1,0),xx);
    table["G0omy0x"]=-G(lst(0,1 - yy,0),xx);
    table["Gomy00x"]=-G(lst(1 - yy,0,0),xx);
    table["Gomyomy0x"]=G(lst(1 - yy,1 - yy,0),xx);
    table["G0omy10x"]=G(lst(0,1 - yy,1,0),xx);
    table["Gomy010x"]=G(lst(1 - yy,0,1,0),xx);
    table["Gomy100x"]=G(lst(1 - yy,1,0,0),xx);
    table["Gomyomy10x"]=-G(lst(1 - yy,1 - yy,1,0),xx);
    table["G0omz0x"]=-G(lst(0,1 - zz,0),xx);
    table["Gomz00x"]=-G(lst(1 - zz,0,0),xx);
    table["Gomzomz0x"]=G(lst(1 - zz,1 - zz,0),xx);
    table["G0omz10x"]=G(lst(0,1 - zz,1,0),xx);
    table["Gomz010x"]=G(lst(1 - zz,0,1,0),xx);
    table["Gomz100x"]=G(lst(1 - zz,1,0,0),xx);
    table["Gomzomz10x"]=-G(lst(1 - zz,1 - zz,1,0),xx);
    table["G0yx"]=G(lst(0,-yy),xx);
    table["G0zx"]=G(lst(0,-zz),xx);
    table["G001x"]=-G(lst(0,0,1),xx);
    table["G011x"]=G(lst(0,1,1),xx);
    table["G101x"]=G(lst(1,0,1),xx);
    table["G111x"]=-G(lst(1,1,1),xx);
    table["G00omyx"]=-G(lst(0,0,1 - yy),xx);
    table["G0omyomyx"]=G(lst(0,1 - yy,1 - yy),xx);
    table["G10omyx"]=G(lst(1,0,1 - yy),xx);
    table["G1omy0x"]=G(lst(1,1 - yy,0),xx);
    table["Gomy0omyx"]=G(lst(1 - yy,0,1 - yy),xx);
    table["Gomyomyomyx"]=-G(lst(1 - yy,1 - yy,1 - yy),xx);
    table["G00omzx"]=-G(lst(0,0,1 - zz),xx);
    table["G0omzomzx"]=G(lst(0,1 - zz,1 - zz),xx);
    table["G10omzx"]=G(lst(1,0,1 - zz),xx);
    table["G1omz0x"]=G(lst(1,1 - zz,0),xx);
    table["Gomz0omzx"]=G(lst(1 - zz,0,1 - zz),xx);
    table["Gomzomzomzx"]=-G(lst(1 - zz,1 - zz,1 - zz),xx);
    table["Gomzzx"]=-G(lst(1 - zz,-zz),xx);
    table["G0zomzx"]=-G(lst(0,-zz,1 - zz),xx);
    table["Gomzzomzx"]=G(lst(1 - zz,-zz,1 - zz),xx);
    table["Gz0x"]=G(lst(-zz,0),xx);
    table["Gzzx"]=G(lst(-zz,-zz),xx);
    table["Gz0omzx"]=-G(lst(-zz,0,1 - zz),xx);
    table["Gzomz0x"]=-G(lst(-zz,1 - zz,0),xx);
    table["Gzomzomzx"]=G(lst(-zz,1 - zz,1 - zz),xx);
    table["Gzzomzx"]=-G(lst(-zz,-zz,1 - zz),xx);
    table["Gomyyx"]=-G(lst(1 - yy,-yy),xx);
    table["G0yomyx"]=-G(lst(0,-yy,1 - yy),xx);
    table["Gomyyomyx"]=G(lst(1 - yy,-yy,1 - yy),xx);
    table["Gy0x"]=G(lst(-yy,0),xx);
    table["Gyyx"]=G(lst(-yy,-yy),xx);
    table["Gy0omyx"]=-G(lst(-yy,0,1 - yy),xx);
    table["Gyomy0x"]=-G(lst(-yy,1 - yy,0),xx);
    table["Gyomyomyx"]=G(lst(-yy,1 - yy,1 - yy),xx);
    table["Gyyomyx"]=-G(lst(-yy,-yy,1 - yy),xx);
    table["G0000x"]=G(lst(0,0,0,0),xx);
    table["Gomy1x"]=G(lst(1 - yy,1),xx);
    table["Gomy110x"]=-G(lst(1 - yy,1,1,0),xx);
    table["Gomz1x"]=G(lst(1 - zz,1),xx);
    table["Gomz110x"]=-G(lst(1 - zz,1,1,0),xx);
    table["G0001x"]=-G(lst(0,0,0,1),xx);
    table["G0101x"]=G(lst(0,1,0,1),xx);
    table["G0111x"]=-G(lst(0,1,1,1),xx);
    table["G0011x"]=G(lst(0,0,1,1),xx);
    table["G1001x"]=G(lst(1,0,0,1),xx);
    table["G1011x"]=-G(lst(1,0,1,1),xx);
    table["G1101x"]=-G(lst(1,1,0,1),xx);
    table["G1111x"]=G(lst(1,1,1,1),xx);
    table["G1omyx"]=G(lst(1,1 - yy),xx);
    table["G1omy10x"]=-G(lst(1,1 - yy,1,0),xx);
    table["Gy010x"]=-G(lst(-yy,0,1,0),xx);
    table["Gyomy10x"]=G(lst(-yy,1 - yy,1,0),xx);
    table["G00yx"]=G(lst(0,0,-yy),xx);
    table["G0omyyx"]=-G(lst(0,1 - yy,-yy),xx);
    table["G0yyx"]=G(lst(0,-yy,-yy),xx);
    table["Gomy0yx"]=-G(lst(1 - yy,0,-yy),xx);
    table["Gomyomyyx"]=G(lst(1 - yy,1 - yy,-yy),xx);
    table["Gomyyyx"]=-G(lst(1 - yy,-yy,-yy),xx);
    table["G00yomyx"]=-G(lst(0,0,-yy,1 - yy),xx);
    table["G010omyx"]=G(lst(0,1,0,1 - yy),xx);
    table["G01omy0x"]=G(lst(0,1,1 - yy,0),xx);
    table["G0omyyomyx"]=G(lst(0,1 - yy,-yy,1 - yy),xx);
    table["G0yomyomyx"]=G(lst(0,-yy,1 - yy,1 - yy),xx);
    table["G0yyomyx"]=-G(lst(0,-yy,-yy,1 - yy),xx);
    table["Gomy0yomyx"]=G(lst(1 - yy,0,-yy,1 - yy),xx);
    table["Gomy10omyx"]=-G(lst(1 - yy,1,0,1 - yy),xx);
    table["Gomy1omy0x"]=-G(lst(1 - yy,1,1 - yy,0),xx);
    table["Gomyomyyomyx"]=-G(lst(1 - yy,1 - yy,-yy,1 - yy),xx);
    table["Gomyyomyomyx"]=-G(lst(1 - yy,-yy,1 - yy,1 - yy),xx);
    table["Gomyyyomyx"]=G(lst(1 - yy,-yy,-yy,1 - yy),xx);
    table["G10yx"]=-G(lst(1,0,-yy),xx);
    table["G1omyyx"]=G(lst(1,1 - yy,-yy),xx);
    table["G000omyx"]=-G(lst(0,0,0,1 - yy),xx);
    table["G00omy0x"]=-G(lst(0,0,1 - yy,0),xx);
    table["G00omyomyx"]=G(lst(0,0,1 - yy,1 - yy),xx);
    table["G0omy00x"]=-G(lst(0,1 - yy,0,0),xx);
    table["G0omy0omyx"]=G(lst(0,1 - yy,0,1 - yy),xx);
    table["G0omyomy0x"]=G(lst(0,1 - yy,1 - yy,0),xx);
    table["G0omyomyomyx"]=-G(lst(0,1 - yy,1 - yy,1 - yy),xx);
    table["G100omyx"]=G(lst(1,0,0,1 - yy),xx);
    table["G10omy0x"]=G(lst(1,0,1 - yy,0),xx);
    table["G10omyomyx"]=-G(lst(1,0,1 - yy,1 - yy),xx);
    table["G10yomyx"]=G(lst(1,0,-yy,1 - yy),xx);
    table["G1omy00x"]=G(lst(1,1 - yy,0,0),xx);
    table["G1omy0omyx"]=-G(lst(1,1 - yy,0,1 - yy),xx);
    table["G1omyomy0x"]=-G(lst(1,1 - yy,1 - yy,0),xx);
    table["G1omyyomyx"]=-G(lst(1,1 - yy,-yy,1 - yy),xx);
    table["Gomy000x"]=-G(lst(1 - yy,0,0,0),xx);
    table["Gomy00omyx"]=G(lst(1 - yy,0,0,1 - yy),xx);
    table["Gomy0omy0x"]=G(lst(1 - yy,0,1 - yy,0),xx);
    table["Gomy0omyomyx"]=-G(lst(1 - yy,0,1 - yy,1 - yy),xx);
    table["Gomyomy00x"]=G(lst(1 - yy,1 - yy,0,0),xx);
    table["Gomyomy0omyx"]=-G(lst(1 - yy,1 - yy,0,1 - yy),xx);
    table["Gomyomyomy0x"]=-G(lst(1 - yy,1 - yy,1 - yy,0),xx);
    table["Gomyomyomyomyx"]=G(lst(1 - yy,1 - yy,1 - yy,1 - yy),xx);
    table["G1omzx"]=G(lst(1,1 - zz),xx);
    table["G1omz10x"]=-G(lst(1,1 - zz,1,0),xx);
    table["Gz010x"]=-G(lst(-zz,0,1,0),xx);
    table["Gzomz10x"]=G(lst(-zz,1 - zz,1,0),xx);
    table["G00zx"]=G(lst(0,0,-zz),xx);
    table["G0omzzx"]=-G(lst(0,1 - zz,-zz),xx);
    table["G0zzx"]=G(lst(0,-zz,-zz),xx);
    table["Gomz0zx"]=-G(lst(1 - zz,0,-zz),xx);
    table["Gomzomzzx"]=G(lst(1 - zz,1 - zz,-zz),xx);
    table["Gomzzzx"]=-G(lst(1 - zz,-zz,-zz),xx);
    table["G00zomzx"]=-G(lst(0,0,-zz,1 - zz),xx);
    table["G010omzx"]=G(lst(0,1,0,1 - zz),xx);
    table["G01omz0x"]=G(lst(0,1,1 - zz,0),xx);
    table["G0omzzomzx"]=G(lst(0,1 - zz,-zz,1 - zz),xx);
    table["G0zomzomzx"]=G(lst(0,-zz,1 - zz,1 - zz),xx);
    table["G0zzomzx"]=-G(lst(0,-zz,-zz,1 - zz),xx);
    table["Gomz0zomzx"]=G(lst(1 - zz,0,-zz,1 - zz),xx);
    table["Gomz10omzx"]=-G(lst(1 - zz,1,0,1 - zz),xx);
    table["Gomz1omz0x"]=-G(lst(1 - zz,1,1 - zz,0),xx);
    table["Gomzomzzomzx"]=-G(lst(1 - zz,1 - zz,-zz,1 - zz),xx);
    table["Gomzzomzomzx"]=-G(lst(1 - zz,-zz,1 - zz,1 - zz),xx);
    table["Gomzzzomzx"]=G(lst(1 - zz,-zz,-zz,1 - zz),xx);
    table["G10zx"]=-G(lst(1,0,-zz),xx);
    table["G1omzzx"]=G(lst(1,1 - zz,-zz),xx);
    table["G000omzx"]=-G(lst(0,0,0,1 - zz),xx);
    table["G00omz0x"]=-G(lst(0,0,1 - zz,0),xx);
    table["G00omzomzx"]=G(lst(0,0,1 - zz,1 - zz),xx);
    table["G0omz00x"]=-G(lst(0,1 - zz,0,0),xx);
    table["G0omz0omzx"]=G(lst(0,1 - zz,0,1 - zz),xx);
    table["G0omzomz0x"]=G(lst(0,1 - zz,1 - zz,0),xx);
    table["G0omzomzomzx"]=-G(lst(0,1 - zz,1 - zz,1 - zz),xx);
    table["G100omzx"]=G(lst(1,0,0,1 - zz),xx);
    table["G10omz0x"]=G(lst(1,0,1 - zz,0),xx);
    table["G10omzomzx"]=-G(lst(1,0,1 - zz,1 - zz),xx);
    table["G10zomzx"]=G(lst(1,0,-zz,1 - zz),xx);
    table["G1omz00x"]=G(lst(1,1 - zz,0,0),xx);
    table["G1omz0omzx"]=-G(lst(1,1 - zz,0,1 - zz),xx);
    table["G1omzomz0x"]=-G(lst(1,1 - zz,1 - zz,0),xx);
    table["G1omzzomzx"]=-G(lst(1,1 - zz,-zz,1 - zz),xx);
    table["Gomz000x"]=-G(lst(1 - zz,0,0,0),xx);
    table["Gomz00omzx"]=G(lst(1 - zz,0,0,1 - zz),xx);
    table["Gomz0omz0x"]=G(lst(1 - zz,0,1 - zz,0),xx);
    table["Gomz0omzomzx"]=-G(lst(1 - zz,0,1 - zz,1 - zz),xx);
    table["Gomzomz00x"]=G(lst(1 - zz,1 - zz,0,0),xx);
    table["Gomzomz0omzx"]=-G(lst(1 - zz,1 - zz,0,1 - zz),xx);
    table["Gomzomzomz0x"]=-G(lst(1 - zz,1 - zz,1 - zz,0),xx);
    table["Gomzomzomzomzx"]=G(lst(1 - zz,1 - zz,1 - zz,1 - zz),xx);
    table["Gomyy0x"]=-G(lst(1 - yy,-yy,0),xx);
    table["Gomyy0omyx"]=G(lst(1 - yy,-yy,0,1 - yy),xx);
    table["Gomyyomy0x"]=G(lst(1 - yy,-yy,1 - yy,0),xx);
    table["Gy0omyomyx"]=G(lst(-yy,0,1 - yy,1 - yy),xx);
    table["Gyomy0omyx"]=G(lst(-yy,1 - yy,0,1 - yy),xx);
    table["Gyomyomy0x"]=G(lst(-yy,1 - yy,1 - yy,0),xx);
    table["G0y0x"]=G(lst(0,-yy,0),xx);
    table["Gy00x"]=G(lst(-yy,0,0),xx);
    table["Gy0yx"]=G(lst(-yy,0,-yy),xx);
    table["Gyyyx"]=G(lst(-yy,-yy,-yy),xx);
    table["G0y0omyx"]=-G(lst(0,-yy,0,1 - yy),xx);
    table["G0yomy0x"]=-G(lst(0,-yy,1 - yy,0),xx);
    table["Gy00omyx"]=-G(lst(-yy,0,0,1 - yy),xx);
    table["Gy0omy0x"]=-G(lst(-yy,0,1 - yy,0),xx);
    table["Gy0yomyx"]=-G(lst(-yy,0,-yy,1 - yy),xx);
    table["Gyomy00x"]=-G(lst(-yy,1 - yy,0,0),xx);
    table["Gyomyomyomyx"]=-G(lst(-yy,1 - yy,1 - yy,1 - yy),xx);
    table["Gyyomyomyx"]=G(lst(-yy,-yy,1 - yy,1 - yy),xx);
    table["Gyyyomyx"]=-G(lst(-yy,-yy,-yy,1 - yy),xx);
    table["Gyomyyx"]=-G(lst(-yy,1 - yy,-yy),xx);
    table["Gyomyyomyx"]=G(lst(-yy,1 - yy,-yy,1 - yy),xx);
    table["Gyy0x"]=G(lst(-yy,-yy,0),xx);
    table["Gyy0omyx"]=-G(lst(-yy,-yy,0,1 - yy),xx);
    table["Gyyomy0x"]=-G(lst(-yy,-yy,1 - yy,0),xx);
    table["Gomzz0x"]=-G(lst(1 - zz,-zz,0),xx);
    table["Gomzz0omzx"]=G(lst(1 - zz,-zz,0,1 - zz),xx);
    table["Gomzzomz0x"]=G(lst(1 - zz,-zz,1 - zz,0),xx);
    table["Gz0omzomzx"]=G(lst(-zz,0,1 - zz,1 - zz),xx);
    table["Gzomz0omzx"]=G(lst(-zz,1 - zz,0,1 - zz),xx);
    table["Gzomzomz0x"]=G(lst(-zz,1 - zz,1 - zz,0),xx);
    table["G0z0x"]=G(lst(0,-zz,0),xx);
    table["Gz00x"]=G(lst(-zz,0,0),xx);
    table["Gz0zx"]=G(lst(-zz,0,-zz),xx);
    table["Gzzzx"]=G(lst(-zz,-zz,-zz),xx);
    table["G0z0omzx"]=-G(lst(0,-zz,0,1 - zz),xx);
    table["G0zomz0x"]=-G(lst(0,-zz,1 - zz,0),xx);
    table["Gz00omzx"]=-G(lst(-zz,0,0,1 - zz),xx);
    table["Gz0omz0x"]=-G(lst(-zz,0,1 - zz,0),xx);
    table["Gz0zomzx"]=-G(lst(-zz,0,-zz,1 - zz),xx);
    table["Gzomz00x"]=-G(lst(-zz,1 - zz,0,0),xx);
    table["Gzomzomzomzx"]=-G(lst(-zz,1 - zz,1 - zz,1 - zz),xx);
    table["Gzzomzomzx"]=G(lst(-zz,-zz,1 - zz,1 - zz),xx);
    table["Gzzzomzx"]=-G(lst(-zz,-zz,-zz,1 - zz),xx);
    table["Gzomzzx"]=-G(lst(-zz,1 - zz,-zz),xx);
    table["Gzomzzomzx"]=G(lst(-zz,1 - zz,-zz,1 - zz),xx);
    table["Gzz0x"]=G(lst(-zz,-zz,0),xx);
    table["Gzz0omzx"]=-G(lst(-zz,-zz,0,1 - zz),xx);
    table["Gzzomz0x"]=-G(lst(-zz,-zz,1 - zz,0),xx);
    table["G0y"]=G(lst(0),yy);
    table["G1y"]=-G(lst(1),yy);
    table["G00y"]=G(lst(0,0),yy);
    table["G10y"]=-G(lst(1,0),yy);
    table["G01y"]=-G(lst(0,1),yy);
    table["G0omxy"]=-G(lst(0,1 - xx),yy);
    table["Gomx0y"]=-G(lst(1 - xx,0),yy);
    table["Gomxomxy"]=G(lst(1 - xx,1 - xx),yy);
    table["Gxomxy"]=-G(lst(-xx,1 - xx),yy);
    table["G11y"]=G(lst(1,1),yy);
    table["G0omzy"]=-G(lst(0,1 - zz),yy);
    table["Gomz0y"]=-G(lst(1 - zz,0),yy);
    table["Gomzomzy"]=G(lst(1 - zz,1 - zz),yy);
    table["Gzomzy"]=-G(lst(-zz,1 - zz),yy);
    table["G000y"]=G(lst(0,0,0),yy);
    table["G010y"]=-G(lst(0,1,0),yy);
    table["Gomx10y"]=G(lst(1 - xx,1,0),yy);
    table["G100y"]=-G(lst(1,0,0),yy);
    table["G110y"]=G(lst(1,1,0),yy);
    table["Gomz10y"]=G(lst(1 - zz,1,0),yy);
    table["G0000y"]=G(lst(0,0,0,0),yy);
    table["G0omx0y"]=-G(lst(0,1 - xx,0),yy);
    table["Gomx00y"]=-G(lst(1 - xx,0,0),yy);
    table["Gomxomx0y"]=G(lst(1 - xx,1 - xx,0),yy);
    table["G0010y"]=-G(lst(0,0,1,0),yy);
    table["G0100y"]=-G(lst(0,1,0,0),yy);
    table["G0omx10y"]=G(lst(0,1 - xx,1,0),yy);
    table["G1000y"]=-G(lst(1,0,0,0),yy);
    table["Gomx010y"]=G(lst(1 - xx,0,1,0),yy);
    table["Gomx100y"]=G(lst(1 - xx,1,0,0),yy);
    table["Gomxomx10y"]=-G(lst(1 - xx,1 - xx,1,0),yy);
    table["G0110y"]=G(lst(0,1,1,0),yy);
    table["G1010y"]=G(lst(1,0,1,0),yy);
    table["G1100y"]=G(lst(1,1,0,0),yy);
    table["G1110y"]=-G(lst(1,1,1,0),yy);
    table["G0omz0y"]=-G(lst(0,1 - zz,0),yy);
    table["Gomz00y"]=-G(lst(1 - zz,0,0),yy);
    table["Gomzomz0y"]=G(lst(1 - zz,1 - zz,0),yy);
    table["G0omz10y"]=G(lst(0,1 - zz,1,0),yy);
    table["Gomz010y"]=G(lst(1 - zz,0,1,0),yy);
    table["Gomz100y"]=G(lst(1 - zz,1,0,0),yy);
    table["Gomzomz10y"]=-G(lst(1 - zz,1 - zz,1,0),yy);
    table["G0xy"]=G(lst(0,-xx),yy);
    table["G0zy"]=G(lst(0,-zz),yy);
    table["G00omxy"]=-G(lst(0,0,1 - xx),yy);
    table["G0omxomxy"]=G(lst(0,1 - xx,1 - xx),yy);
    table["G10omxy"]=G(lst(1,0,1 - xx),yy);
    table["G1omx0y"]=G(lst(1,1 - xx,0),yy);
    table["Gomx0omxy"]=G(lst(1 - xx,0,1 - xx),yy);
    table["Gomxomxomxy"]=-G(lst(1 - xx,1 - xx,1 - xx),yy);
    table["G001y"]=-G(lst(0,0,1),yy);
    table["G011y"]=G(lst(0,1,1),yy);
    table["G101y"]=G(lst(1,0,1),yy);
    table["G111y"]=-G(lst(1,1,1),yy);
    table["G00omzy"]=-G(lst(0,0,1 - zz),yy);
    table["G0omzomzy"]=G(lst(0,1 - zz,1 - zz),yy);
    table["G10omzy"]=G(lst(1,0,1 - zz),yy);
    table["G1omz0y"]=G(lst(1,1 - zz,0),yy);
    table["Gomz0omzy"]=G(lst(1 - zz,0,1 - zz),yy);
    table["Gomzomzomzy"]=-G(lst(1 - zz,1 - zz,1 - zz),yy);
    table["Gomxxy"]=-G(lst(1 - xx,-xx),yy);
    table["Gx0y"]=G(lst(-xx,0),yy);
    table["Gxxy"]=G(lst(-xx,-xx),yy);
    table["G0xomxy"]=-G(lst(0,-xx,1 - xx),yy);
    table["Gomxxomxy"]=G(lst(1 - xx,-xx,1 - xx),yy);
    table["Gx0omxy"]=-G(lst(-xx,0,1 - xx),yy);
    table["Gxomx0y"]=-G(lst(-xx,1 - xx,0),yy);
    table["Gxomxomxy"]=G(lst(-xx,1 - xx,1 - xx),yy);
    table["Gxxomxy"]=-G(lst(-xx,-xx,1 - xx),yy);
    table["Gomzzy"]=-G(lst(1 - zz,-zz),yy);
    table["G0zomzy"]=-G(lst(0,-zz,1 - zz),yy);
    table["Gomzzomzy"]=G(lst(1 - zz,-zz,1 - zz),yy);
    table["Gomx1y"]=G(lst(1 - xx,1),yy);
    table["Gomx110y"]=-G(lst(1 - xx,1,1,0),yy);
    table["Gomz1y"]=G(lst(1 - zz,1),yy);
    table["Gomz110y"]=-G(lst(1 - zz,1,1,0),yy);
    table["Gx010y"]=-G(lst(-xx,0,1,0),yy);
    table["Gxomx10y"]=G(lst(-xx,1 - xx,1,0),yy);
    table["G00xy"]=G(lst(0,0,-xx),yy);
    table["G0omxxy"]=-G(lst(0,1 - xx,-xx),yy);
    table["G0xxy"]=G(lst(0,-xx,-xx),yy);
    table["Gomx0xy"]=-G(lst(1 - xx,0,-xx),yy);
    table["Gomxomxxy"]=G(lst(1 - xx,1 - xx,-xx),yy);
    table["Gomxxxy"]=-G(lst(1 - xx,-xx,-xx),yy);
    table["G00xomxy"]=-G(lst(0,0,-xx,1 - xx),yy);
    table["G010omxy"]=G(lst(0,1,0,1 - xx),yy);
    table["G01omx0y"]=G(lst(0,1,1 - xx,0),yy);
    table["G0omxxomxy"]=G(lst(0,1 - xx,-xx,1 - xx),yy);
    table["G0xomxomxy"]=G(lst(0,-xx,1 - xx,1 - xx),yy);
    table["G0xxomxy"]=-G(lst(0,-xx,-xx,1 - xx),yy);
    table["Gomx0xomxy"]=G(lst(1 - xx,0,-xx,1 - xx),yy);
    table["Gomx10omxy"]=-G(lst(1 - xx,1,0,1 - xx),yy);
    table["Gomx1omx0y"]=-G(lst(1 - xx,1,1 - xx,0),yy);
    table["Gomxomxxomxy"]=-G(lst(1 - xx,1 - xx,-xx,1 - xx),yy);
    table["Gomxxomxomxy"]=-G(lst(1 - xx,-xx,1 - xx,1 - xx),yy);
    table["Gomxxxomxy"]=G(lst(1 - xx,-xx,-xx,1 - xx),yy);
    table["G1omxy"]=G(lst(1,1 - xx),yy);
    table["G10xy"]=-G(lst(1,0,-xx),yy);
    table["G1omxxy"]=G(lst(1,1 - xx,-xx),yy);
    table["G000omxy"]=-G(lst(0,0,0,1 - xx),yy);
    table["G00omx0y"]=-G(lst(0,0,1 - xx,0),yy);
    table["G00omxomxy"]=G(lst(0,0,1 - xx,1 - xx),yy);
    table["G0omx00y"]=-G(lst(0,1 - xx,0,0),yy);
    table["G0omx0omxy"]=G(lst(0,1 - xx,0,1 - xx),yy);
    table["G0omxomx0y"]=G(lst(0,1 - xx,1 - xx,0),yy);
    table["G0omxomxomxy"]=-G(lst(0,1 - xx,1 - xx,1 - xx),yy);
    table["G100omxy"]=G(lst(1,0,0,1 - xx),yy);
    table["G10omx0y"]=G(lst(1,0,1 - xx,0),yy);
    table["G10omxomxy"]=-G(lst(1,0,1 - xx,1 - xx),yy);
    table["G10xomxy"]=G(lst(1,0,-xx,1 - xx),yy);
    table["G1omx00y"]=G(lst(1,1 - xx,0,0),yy);
    table["G1omx0omxy"]=-G(lst(1,1 - xx,0,1 - xx),yy);
    table["G1omx10y"]=-G(lst(1,1 - xx,1,0),yy);
    table["G1omxomx0y"]=-G(lst(1,1 - xx,1 - xx,0),yy);
    table["G1omxxomxy"]=-G(lst(1,1 - xx,-xx,1 - xx),yy);
    table["Gomx000y"]=-G(lst(1 - xx,0,0,0),yy);
    table["Gomx00omxy"]=G(lst(1 - xx,0,0,1 - xx),yy);
    table["Gomx0omx0y"]=G(lst(1 - xx,0,1 - xx,0),yy);
    table["Gomx0omxomxy"]=-G(lst(1 - xx,0,1 - xx,1 - xx),yy);
    table["Gomxomx00y"]=G(lst(1 - xx,1 - xx,0,0),yy);
    table["Gomxomx0omxy"]=-G(lst(1 - xx,1 - xx,0,1 - xx),yy);
    table["Gomxomxomx0y"]=-G(lst(1 - xx,1 - xx,1 - xx,0),yy);
    table["Gomxomxomxomxy"]=G(lst(1 - xx,1 - xx,1 - xx,1 - xx),yy);
    table["G0001y"]=-G(lst(0,0,0,1),yy);
    table["G0101y"]=G(lst(0,1,0,1),yy);
    table["G0111y"]=-G(lst(0,1,1,1),yy);
    table["G0011y"]=G(lst(0,0,1,1),yy);
    table["G1001y"]=G(lst(1,0,0,1),yy);
    table["G1011y"]=-G(lst(1,0,1,1),yy);
    table["G1101y"]=-G(lst(1,1,0,1),yy);
    table["G1111y"]=G(lst(1,1,1,1),yy);
    table["G1omzy"]=G(lst(1,1 - zz),yy);
    table["G1omz10y"]=-G(lst(1,1 - zz,1,0),yy);
    table["Gz0y"]=G(lst(-zz,0),yy);
    table["Gzomz0y"]=-G(lst(-zz,1 - zz,0),yy);
    table["Gz010y"]=-G(lst(-zz,0,1,0),yy);
    table["Gzomz10y"]=G(lst(-zz,1 - zz,1,0),yy);
    table["G00zy"]=G(lst(0,0,-zz),yy);
    table["G0omzzy"]=-G(lst(0,1 - zz,-zz),yy);
    table["G0zzy"]=G(lst(0,-zz,-zz),yy);
    table["Gomz0zy"]=-G(lst(1 - zz,0,-zz),yy);
    table["Gomzomzzy"]=G(lst(1 - zz,1 - zz,-zz),yy);
    table["Gomzzzy"]=-G(lst(1 - zz,-zz,-zz),yy);
    table["G00zomzy"]=-G(lst(0,0,-zz,1 - zz),yy);
    table["G010omzy"]=G(lst(0,1,0,1 - zz),yy);
    table["G01omz0y"]=G(lst(0,1,1 - zz,0),yy);
    table["G0omzzomzy"]=G(lst(0,1 - zz,-zz,1 - zz),yy);
    table["G0zomzomzy"]=G(lst(0,-zz,1 - zz,1 - zz),yy);
    table["G0zzomzy"]=-G(lst(0,-zz,-zz,1 - zz),yy);
    table["Gomz0zomzy"]=G(lst(1 - zz,0,-zz,1 - zz),yy);
    table["Gomz10omzy"]=-G(lst(1 - zz,1,0,1 - zz),yy);
    table["Gomz1omz0y"]=-G(lst(1 - zz,1,1 - zz,0),yy);
    table["Gomzomzzomzy"]=-G(lst(1 - zz,1 - zz,-zz,1 - zz),yy);
    table["Gomzzomzomzy"]=-G(lst(1 - zz,-zz,1 - zz,1 - zz),yy);
    table["Gomzzzomzy"]=G(lst(1 - zz,-zz,-zz,1 - zz),yy);
    table["G10zy"]=-G(lst(1,0,-zz),yy);
    table["G1omzzy"]=G(lst(1,1 - zz,-zz),yy);
    table["G000omzy"]=-G(lst(0,0,0,1 - zz),yy);
    table["G00omz0y"]=-G(lst(0,0,1 - zz,0),yy);
    table["G00omzomzy"]=G(lst(0,0,1 - zz,1 - zz),yy);
    table["G0omz00y"]=-G(lst(0,1 - zz,0,0),yy);
    table["G0omz0omzy"]=G(lst(0,1 - zz,0,1 - zz),yy);
    table["G0omzomz0y"]=G(lst(0,1 - zz,1 - zz,0),yy);
    table["G0omzomzomzy"]=-G(lst(0,1 - zz,1 - zz,1 - zz),yy);
    table["G100omzy"]=G(lst(1,0,0,1 - zz),yy);
    table["G10omz0y"]=G(lst(1,0,1 - zz,0),yy);
    table["G10omzomzy"]=-G(lst(1,0,1 - zz,1 - zz),yy);
    table["G10zomzy"]=G(lst(1,0,-zz,1 - zz),yy);
    table["G1omz00y"]=G(lst(1,1 - zz,0,0),yy);
    table["G1omz0omzy"]=-G(lst(1,1 - zz,0,1 - zz),yy);
    table["G1omzomz0y"]=-G(lst(1,1 - zz,1 - zz,0),yy);
    table["G1omzzomzy"]=-G(lst(1,1 - zz,-zz,1 - zz),yy);
    table["Gomz000y"]=-G(lst(1 - zz,0,0,0),yy);
    table["Gomz00omzy"]=G(lst(1 - zz,0,0,1 - zz),yy);
    table["Gomz0omz0y"]=G(lst(1 - zz,0,1 - zz,0),yy);
    table["Gomz0omzomzy"]=-G(lst(1 - zz,0,1 - zz,1 - zz),yy);
    table["Gomzomz00y"]=G(lst(1 - zz,1 - zz,0,0),yy);
    table["Gomzomz0omzy"]=-G(lst(1 - zz,1 - zz,0,1 - zz),yy);
    table["Gomzomzomz0y"]=-G(lst(1 - zz,1 - zz,1 - zz,0),yy);
    table["Gomzomzomzomzy"]=G(lst(1 - zz,1 - zz,1 - zz,1 - zz),yy);
    table["Gzzy"]=G(lst(-zz,-zz),yy);
    table["Gz0omzy"]=-G(lst(-zz,0,1 - zz),yy);
    table["Gzomzomzy"]=G(lst(-zz,1 - zz,1 - zz),yy);
    table["Gzzomzy"]=-G(lst(-zz,-zz,1 - zz),yy);
    table["Gomxx0y"]=-G(lst(1 - xx,-xx,0),yy);
    table["Gomxx0omxy"]=G(lst(1 - xx,-xx,0,1 - xx),yy);
    table["Gomxxomx0y"]=G(lst(1 - xx,-xx,1 - xx,0),yy);
    table["Gx0omxomxy"]=G(lst(-xx,0,1 - xx,1 - xx),yy);
    table["Gxomx0omxy"]=G(lst(-xx,1 - xx,0,1 - xx),yy);
    table["Gxomxomx0y"]=G(lst(-xx,1 - xx,1 - xx,0),yy);
    table["G0x0y"]=G(lst(0,-xx,0),yy);
    table["Gx00y"]=G(lst(-xx,0,0),yy);
    table["Gx0xy"]=G(lst(-xx,0,-xx),yy);
    table["Gxxxy"]=G(lst(-xx,-xx,-xx),yy);
    table["G0x0omxy"]=-G(lst(0,-xx,0,1 - xx),yy);
    table["G0xomx0y"]=-G(lst(0,-xx,1 - xx,0),yy);
    table["Gx00omxy"]=-G(lst(-xx,0,0,1 - xx),yy);
    table["Gx0omx0y"]=-G(lst(-xx,0,1 - xx,0),yy);
    table["Gx0xomxy"]=-G(lst(-xx,0,-xx,1 - xx),yy);
    table["Gxomx00y"]=-G(lst(-xx,1 - xx,0,0),yy);
    table["Gxomxomxomxy"]=-G(lst(-xx,1 - xx,1 - xx,1 - xx),yy);
    table["Gxxomxomxy"]=G(lst(-xx,-xx,1 - xx,1 - xx),yy);
    table["Gxxxomxy"]=-G(lst(-xx,-xx,-xx,1 - xx),yy);
    table["Gxomxxy"]=-G(lst(-xx,1 - xx,-xx),yy);
    table["Gxomxxomxy"]=G(lst(-xx,1 - xx,-xx,1 - xx),yy);
    table["Gxx0y"]=G(lst(-xx,-xx,0),yy);
    table["Gxx0omxy"]=-G(lst(-xx,-xx,0,1 - xx),yy);
    table["Gxxomx0y"]=-G(lst(-xx,-xx,1 - xx,0),yy);
    table["Gomzz0y"]=-G(lst(1 - zz,-zz,0),yy);
    table["Gomzz0omzy"]=G(lst(1 - zz,-zz,0,1 - zz),yy);
    table["Gomzzomz0y"]=G(lst(1 - zz,-zz,1 - zz,0),yy);
    table["Gz0omzomzy"]=G(lst(-zz,0,1 - zz,1 - zz),yy);
    table["Gzomz0omzy"]=G(lst(-zz,1 - zz,0,1 - zz),yy);
    table["Gzomzomz0y"]=G(lst(-zz,1 - zz,1 - zz,0),yy);
    table["G0z0y"]=G(lst(0,-zz,0),yy);
    table["Gz00y"]=G(lst(-zz,0,0),yy);
    table["Gz0zy"]=G(lst(-zz,0,-zz),yy);
    table["Gzzzy"]=G(lst(-zz,-zz,-zz),yy);
    table["G0z0omzy"]=-G(lst(0,-zz,0,1 - zz),yy);
    table["G0zomz0y"]=-G(lst(0,-zz,1 - zz,0),yy);
    table["Gz00omzy"]=-G(lst(-zz,0,0,1 - zz),yy);
    table["Gz0omz0y"]=-G(lst(-zz,0,1 - zz,0),yy);
    table["Gz0zomzy"]=-G(lst(-zz,0,-zz,1 - zz),yy);
    table["Gzomz00y"]=-G(lst(-zz,1 - zz,0,0),yy);
    table["Gzomzomzomzy"]=-G(lst(-zz,1 - zz,1 - zz,1 - zz),yy);
    table["Gzzomzomzy"]=G(lst(-zz,-zz,1 - zz,1 - zz),yy);
    table["Gzzzomzy"]=-G(lst(-zz,-zz,-zz,1 - zz),yy);
    table["Gzomzzy"]=-G(lst(-zz,1 - zz,-zz),yy);
    table["Gzomzzomzy"]=G(lst(-zz,1 - zz,-zz,1 - zz),yy);
    table["Gzz0y"]=G(lst(-zz,-zz,0),yy);
    table["Gzz0omzy"]=-G(lst(-zz,-zz,0,1 - zz),yy);
    table["Gzzomz0y"]=-G(lst(-zz,-zz,1 - zz,0),yy);
    table["G0z"]=G(lst(0),zz);
    table["G1z"]=-G(lst(1),zz);
    table["G00z"]=G(lst(0,0),zz);
    table["G10z"]=-G(lst(1,0),zz);
    table["G01z"]=-G(lst(0,1),zz);
    table["G0omxz"]=-G(lst(0,1 - xx),zz);
    table["Gomx0z"]=-G(lst(1 - xx,0),zz);
    table["Gomxomxz"]=G(lst(1 - xx,1 - xx),zz);
    table["Gxomxz"]=-G(lst(-xx,1 - xx),zz);
    table["G0omyz"]=-G(lst(0,1 - yy),zz);
    table["Gomy0z"]=-G(lst(1 - yy,0),zz);
    table["Gomyomyz"]=G(lst(1 - yy,1 - yy),zz);
    table["Gyomyz"]=-G(lst(-yy,1 - yy),zz);
    table["G11z"]=G(lst(1,1),zz);
    table["G000z"]=G(lst(0,0,0),zz);
    table["G010z"]=-G(lst(0,1,0),zz);
    table["Gomx10z"]=G(lst(1 - xx,1,0),zz);
    table["G100z"]=-G(lst(1,0,0),zz);
    table["Gomy10z"]=G(lst(1 - yy,1,0),zz);
    table["G110z"]=G(lst(1,1,0),zz);
    table["G0000z"]=G(lst(0,0,0,0),zz);
    table["G0omx0z"]=-G(lst(0,1 - xx,0),zz);
    table["Gomx00z"]=-G(lst(1 - xx,0,0),zz);
    table["Gomxomx0z"]=G(lst(1 - xx,1 - xx,0),zz);
    table["G0010z"]=-G(lst(0,0,1,0),zz);
    table["G0100z"]=-G(lst(0,1,0,0),zz);
    table["G0omx10z"]=G(lst(0,1 - xx,1,0),zz);
    table["G1000z"]=-G(lst(1,0,0,0),zz);
    table["Gomx010z"]=G(lst(1 - xx,0,1,0),zz);
    table["Gomx100z"]=G(lst(1 - xx,1,0,0),zz);
    table["Gomxomx10z"]=-G(lst(1 - xx,1 - xx,1,0),zz);
    table["G0omy0z"]=-G(lst(0,1 - yy,0),zz);
    table["Gomy00z"]=-G(lst(1 - yy,0,0),zz);
    table["Gomyomy0z"]=G(lst(1 - yy,1 - yy,0),zz);
    table["G0omy10z"]=G(lst(0,1 - yy,1,0),zz);
    table["Gomy010z"]=G(lst(1 - yy,0,1,0),zz);
    table["Gomy100z"]=G(lst(1 - yy,1,0,0),zz);
    table["Gomyomy10z"]=-G(lst(1 - yy,1 - yy,1,0),zz);
    table["G0110z"]=G(lst(0,1,1,0),zz);
    table["G1010z"]=G(lst(1,0,1,0),zz);
    table["G1100z"]=G(lst(1,1,0,0),zz);
    table["G1110z"]=-G(lst(1,1,1,0),zz);
    table["G0xz"]=G(lst(0,-xx),zz);
    table["G0yz"]=G(lst(0,-yy),zz);
    table["G00omxz"]=-G(lst(0,0,1 - xx),zz);
    table["G0omxomxz"]=G(lst(0,1 - xx,1 - xx),zz);
    table["G10omxz"]=G(lst(1,0,1 - xx),zz);
    table["G1omx0z"]=G(lst(1,1 - xx,0),zz);
    table["Gomx0omxz"]=G(lst(1 - xx,0,1 - xx),zz);
    table["Gomxomxomxz"]=-G(lst(1 - xx,1 - xx,1 - xx),zz);
    table["G00omyz"]=-G(lst(0,0,1 - yy),zz);
    table["G0omyomyz"]=G(lst(0,1 - yy,1 - yy),zz);
    table["G10omyz"]=G(lst(1,0,1 - yy),zz);
    table["G1omy0z"]=G(lst(1,1 - yy,0),zz);
    table["Gomy0omyz"]=G(lst(1 - yy,0,1 - yy),zz);
    table["Gomyomyomyz"]=-G(lst(1 - yy,1 - yy,1 - yy),zz);
    table["G001z"]=-G(lst(0,0,1),zz);
    table["G011z"]=G(lst(0,1,1),zz);
    table["G101z"]=G(lst(1,0,1),zz);
    table["G111z"]=-G(lst(1,1,1),zz);
    table["Gomxxz"]=-G(lst(1 - xx,-xx),zz);
    table["Gx0z"]=G(lst(-xx,0),zz);
    table["Gxxz"]=G(lst(-xx,-xx),zz);
    table["G0xomxz"]=-G(lst(0,-xx,1 - xx),zz);
    table["Gomxxomxz"]=G(lst(1 - xx,-xx,1 - xx),zz);
    table["Gx0omxz"]=-G(lst(-xx,0,1 - xx),zz);
    table["Gxomx0z"]=-G(lst(-xx,1 - xx,0),zz);
    table["Gxomxomxz"]=G(lst(-xx,1 - xx,1 - xx),zz);
    table["Gxxomxz"]=-G(lst(-xx,-xx,1 - xx),zz);
    table["Gomyyz"]=-G(lst(1 - yy,-yy),zz);
    table["G0yomyz"]=-G(lst(0,-yy,1 - yy),zz);
    table["Gomyyomyz"]=G(lst(1 - yy,-yy,1 - yy),zz);
    table["Gy0z"]=G(lst(-yy,0),zz);
    table["Gyyz"]=G(lst(-yy,-yy),zz);
    table["Gy0omyz"]=-G(lst(-yy,0,1 - yy),zz);
    table["Gyomy0z"]=-G(lst(-yy,1 - yy,0),zz);
    table["Gyomyomyz"]=G(lst(-yy,1 - yy,1 - yy),zz);
    table["Gyyomyz"]=-G(lst(-yy,-yy,1 - yy),zz);
    table["Gomx1z"]=G(lst(1 - xx,1),zz);
    table["Gomx110z"]=-G(lst(1 - xx,1,1,0),zz);
    table["Gomy1z"]=G(lst(1 - yy,1),zz);
    table["Gomy110z"]=-G(lst(1 - yy,1,1,0),zz);
    table["Gx010z"]=-G(lst(-xx,0,1,0),zz);
    table["Gxomx10z"]=G(lst(-xx,1 - xx,1,0),zz);
    table["G00xz"]=G(lst(0,0,-xx),zz);
    table["G0omxxz"]=-G(lst(0,1 - xx,-xx),zz);
    table["G0xxz"]=G(lst(0,-xx,-xx),zz);
    table["Gomx0xz"]=-G(lst(1 - xx,0,-xx),zz);
    table["Gomxomxxz"]=G(lst(1 - xx,1 - xx,-xx),zz);
    table["Gomxxxz"]=-G(lst(1 - xx,-xx,-xx),zz);
    table["G00xomxz"]=-G(lst(0,0,-xx,1 - xx),zz);
    table["G010omxz"]=G(lst(0,1,0,1 - xx),zz);
    table["G01omx0z"]=G(lst(0,1,1 - xx,0),zz);
    table["G0omxxomxz"]=G(lst(0,1 - xx,-xx,1 - xx),zz);
    table["G0xomxomxz"]=G(lst(0,-xx,1 - xx,1 - xx),zz);
    table["G0xxomxz"]=-G(lst(0,-xx,-xx,1 - xx),zz);
    table["Gomx0xomxz"]=G(lst(1 - xx,0,-xx,1 - xx),zz);
    table["Gomx10omxz"]=-G(lst(1 - xx,1,0,1 - xx),zz);
    table["Gomx1omx0z"]=-G(lst(1 - xx,1,1 - xx,0),zz);
    table["Gomxomxxomxz"]=-G(lst(1 - xx,1 - xx,-xx,1 - xx),zz);
    table["Gomxxomxomxz"]=-G(lst(1 - xx,-xx,1 - xx,1 - xx),zz);
    table["Gomxxxomxz"]=G(lst(1 - xx,-xx,-xx,1 - xx),zz);
    table["Gy010z"]=-G(lst(-yy,0,1,0),zz);
    table["Gyomy10z"]=G(lst(-yy,1 - yy,1,0),zz);
    table["G00yz"]=G(lst(0,0,-yy),zz);
    table["G0omyyz"]=-G(lst(0,1 - yy,-yy),zz);
    table["G0yyz"]=G(lst(0,-yy,-yy),zz);
    table["Gomy0yz"]=-G(lst(1 - yy,0,-yy),zz);
    table["Gomyomyyz"]=G(lst(1 - yy,1 - yy,-yy),zz);
    table["Gomyyyz"]=-G(lst(1 - yy,-yy,-yy),zz);
    table["G00yomyz"]=-G(lst(0,0,-yy,1 - yy),zz);
    table["G010omyz"]=G(lst(0,1,0,1 - yy),zz);
    table["G01omy0z"]=G(lst(0,1,1 - yy,0),zz);
    table["G0omyyomyz"]=G(lst(0,1 - yy,-yy,1 - yy),zz);
    table["G0yomyomyz"]=G(lst(0,-yy,1 - yy,1 - yy),zz);
    table["G0yyomyz"]=-G(lst(0,-yy,-yy,1 - yy),zz);
    table["Gomy0yomyz"]=G(lst(1 - yy,0,-yy,1 - yy),zz);
    table["Gomy10omyz"]=-G(lst(1 - yy,1,0,1 - yy),zz);
    table["Gomy1omy0z"]=-G(lst(1 - yy,1,1 - yy,0),zz);
    table["Gomyomyyomyz"]=-G(lst(1 - yy,1 - yy,-yy,1 - yy),zz);
    table["Gomyyomyomyz"]=-G(lst(1 - yy,-yy,1 - yy,1 - yy),zz);
    table["Gomyyyomyz"]=G(lst(1 - yy,-yy,-yy,1 - yy),zz);
    table["G1omxz"]=G(lst(1,1 - xx),zz);
    table["G10xz"]=-G(lst(1,0,-xx),zz);
    table["G1omxxz"]=G(lst(1,1 - xx,-xx),zz);
    table["G000omxz"]=-G(lst(0,0,0,1 - xx),zz);
    table["G00omx0z"]=-G(lst(0,0,1 - xx,0),zz);
    table["G00omxomxz"]=G(lst(0,0,1 - xx,1 - xx),zz);
    table["G0omx00z"]=-G(lst(0,1 - xx,0,0),zz);
    table["G0omx0omxz"]=G(lst(0,1 - xx,0,1 - xx),zz);
    table["G0omxomx0z"]=G(lst(0,1 - xx,1 - xx,0),zz);
    table["G0omxomxomxz"]=-G(lst(0,1 - xx,1 - xx,1 - xx),zz);
    table["G100omxz"]=G(lst(1,0,0,1 - xx),zz);
    table["G10omx0z"]=G(lst(1,0,1 - xx,0),zz);
    table["G10omxomxz"]=-G(lst(1,0,1 - xx,1 - xx),zz);
    table["G10xomxz"]=G(lst(1,0,-xx,1 - xx),zz);
    table["G1omx00z"]=G(lst(1,1 - xx,0,0),zz);
    table["G1omx0omxz"]=-G(lst(1,1 - xx,0,1 - xx),zz);
    table["G1omx10z"]=-G(lst(1,1 - xx,1,0),zz);
    table["G1omxomx0z"]=-G(lst(1,1 - xx,1 - xx,0),zz);
    table["G1omxxomxz"]=-G(lst(1,1 - xx,-xx,1 - xx),zz);
    table["Gomx000z"]=-G(lst(1 - xx,0,0,0),zz);
    table["Gomx00omxz"]=G(lst(1 - xx,0,0,1 - xx),zz);
    table["Gomx0omx0z"]=G(lst(1 - xx,0,1 - xx,0),zz);
    table["Gomx0omxomxz"]=-G(lst(1 - xx,0,1 - xx,1 - xx),zz);
    table["Gomxomx00z"]=G(lst(1 - xx,1 - xx,0,0),zz);
    table["Gomxomx0omxz"]=-G(lst(1 - xx,1 - xx,0,1 - xx),zz);
    table["Gomxomxomx0z"]=-G(lst(1 - xx,1 - xx,1 - xx,0),zz);
    table["Gomxomxomxomxz"]=G(lst(1 - xx,1 - xx,1 - xx,1 - xx),zz);
    table["G1omyz"]=G(lst(1,1 - yy),zz);
    table["G10yz"]=-G(lst(1,0,-yy),zz);
    table["G1omyyz"]=G(lst(1,1 - yy,-yy),zz);
    table["G000omyz"]=-G(lst(0,0,0,1 - yy),zz);
    table["G00omy0z"]=-G(lst(0,0,1 - yy,0),zz);
    table["G00omyomyz"]=G(lst(0,0,1 - yy,1 - yy),zz);
    table["G0omy00z"]=-G(lst(0,1 - yy,0,0),zz);
    table["G0omy0omyz"]=G(lst(0,1 - yy,0,1 - yy),zz);
    table["G0omyomy0z"]=G(lst(0,1 - yy,1 - yy,0),zz);
    table["G0omyomyomyz"]=-G(lst(0,1 - yy,1 - yy,1 - yy),zz);
    table["G100omyz"]=G(lst(1,0,0,1 - yy),zz);
    table["G10omy0z"]=G(lst(1,0,1 - yy,0),zz);
    table["G10omyomyz"]=-G(lst(1,0,1 - yy,1 - yy),zz);
    table["G10yomyz"]=G(lst(1,0,-yy,1 - yy),zz);
    table["G1omy00z"]=G(lst(1,1 - yy,0,0),zz);
    table["G1omy0omyz"]=-G(lst(1,1 - yy,0,1 - yy),zz);
    table["G1omy10z"]=-G(lst(1,1 - yy,1,0),zz);
    table["G1omyomy0z"]=-G(lst(1,1 - yy,1 - yy,0),zz);
    table["G1omyyomyz"]=-G(lst(1,1 - yy,-yy,1 - yy),zz);
    table["Gomy000z"]=-G(lst(1 - yy,0,0,0),zz);
    table["Gomy00omyz"]=G(lst(1 - yy,0,0,1 - yy),zz);
    table["Gomy0omy0z"]=G(lst(1 - yy,0,1 - yy,0),zz);
    table["Gomy0omyomyz"]=-G(lst(1 - yy,0,1 - yy,1 - yy),zz);
    table["Gomyomy00z"]=G(lst(1 - yy,1 - yy,0,0),zz);
    table["Gomyomy0omyz"]=-G(lst(1 - yy,1 - yy,0,1 - yy),zz);
    table["Gomyomyomy0z"]=-G(lst(1 - yy,1 - yy,1 - yy,0),zz);
    table["Gomyomyomyomyz"]=G(lst(1 - yy,1 - yy,1 - yy,1 - yy),zz);
    table["G0001z"]=-G(lst(0,0,0,1),zz);
    table["G0101z"]=G(lst(0,1,0,1),zz);
    table["G0111z"]=-G(lst(0,1,1,1),zz);
    table["G0011z"]=G(lst(0,0,1,1),zz);
    table["G1001z"]=G(lst(1,0,0,1),zz);
    table["G1011z"]=-G(lst(1,0,1,1),zz);
    table["G1101z"]=-G(lst(1,1,0,1),zz);
    table["G1111z"]=G(lst(1,1,1,1),zz);
    table["Gomxx0z"]=-G(lst(1 - xx,-xx,0),zz);
    table["Gomxx0omxz"]=G(lst(1 - xx,-xx,0,1 - xx),zz);
    table["Gomxxomx0z"]=G(lst(1 - xx,-xx,1 - xx,0),zz);
    table["Gx0omxomxz"]=G(lst(-xx,0,1 - xx,1 - xx),zz);
    table["Gxomx0omxz"]=G(lst(-xx,1 - xx,0,1 - xx),zz);
    table["Gxomxomx0z"]=G(lst(-xx,1 - xx,1 - xx,0),zz);
    table["G0x0z"]=G(lst(0,-xx,0),zz);
    table["Gx00z"]=G(lst(-xx,0,0),zz);
    table["Gx0xz"]=G(lst(-xx,0,-xx),zz);
    table["Gxxxz"]=G(lst(-xx,-xx,-xx),zz);
    table["G0x0omxz"]=-G(lst(0,-xx,0,1 - xx),zz);
    table["G0xomx0z"]=-G(lst(0,-xx,1 - xx,0),zz);
    table["Gx00omxz"]=-G(lst(-xx,0,0,1 - xx),zz);
    table["Gx0omx0z"]=-G(lst(-xx,0,1 - xx,0),zz);
    table["Gx0xomxz"]=-G(lst(-xx,0,-xx,1 - xx),zz);
    table["Gxomx00z"]=-G(lst(-xx,1 - xx,0,0),zz);
    table["Gxomxomxomxz"]=-G(lst(-xx,1 - xx,1 - xx,1 - xx),zz);
    table["Gxxomxomxz"]=G(lst(-xx,-xx,1 - xx,1 - xx),zz);
    table["Gxxxomxz"]=-G(lst(-xx,-xx,-xx,1 - xx),zz);
    table["Gxomxxz"]=-G(lst(-xx,1 - xx,-xx),zz);
    table["Gxomxxomxz"]=G(lst(-xx,1 - xx,-xx,1 - xx),zz);
    table["Gomyy0z"]=-G(lst(1 - yy,-yy,0),zz);
    table["Gomyy0omyz"]=G(lst(1 - yy,-yy,0,1 - yy),zz);
    table["Gomyyomy0z"]=G(lst(1 - yy,-yy,1 - yy,0),zz);
    table["Gy0omyomyz"]=G(lst(-yy,0,1 - yy,1 - yy),zz);
    table["Gyomy0omyz"]=G(lst(-yy,1 - yy,0,1 - yy),zz);
    table["Gyomyomy0z"]=G(lst(-yy,1 - yy,1 - yy,0),zz);
    table["G0y0z"]=G(lst(0,-yy,0),zz);
    table["Gy00z"]=G(lst(-yy,0,0),zz);
    table["Gy0yz"]=G(lst(-yy,0,-yy),zz);
    table["Gyyyz"]=G(lst(-yy,-yy,-yy),zz);
    table["G0y0omyz"]=-G(lst(0,-yy,0,1 - yy),zz);
    table["G0yomy0z"]=-G(lst(0,-yy,1 - yy,0),zz);
    table["Gy00omyz"]=-G(lst(-yy,0,0,1 - yy),zz);
    table["Gy0omy0z"]=-G(lst(-yy,0,1 - yy,0),zz);
    table["Gy0yomyz"]=-G(lst(-yy,0,-yy,1 - yy),zz);
    table["Gyomy00z"]=-G(lst(-yy,1 - yy,0,0),zz);
    table["Gyomyomyomyz"]=-G(lst(-yy,1 - yy,1 - yy,1 - yy),zz);
    table["Gyyomyomyz"]=G(lst(-yy,-yy,1 - yy,1 - yy),zz);
    table["Gyyyomyz"]=-G(lst(-yy,-yy,-yy,1 - yy),zz);
    table["Gxx0z"]=G(lst(-xx,-xx,0),zz);
    table["Gxx0omxz"]=-G(lst(-xx,-xx,0,1 - xx),zz);
    table["Gxxomx0z"]=-G(lst(-xx,-xx,1 - xx,0),zz);
    table["Gyomyyz"]=-G(lst(-yy,1 - yy,-yy),zz);
    table["Gyy0z"]=G(lst(-yy,-yy,0),zz);
    table["Gyomyyomyz"]=G(lst(-yy,1 - yy,-yy,1 - yy),zz);
    table["Gyy0omyz"]=-G(lst(-yy,-yy,0,1 - yy),zz);
    table["Gyyomy0z"]=-G(lst(-yy,-yy,1 - yy,0),zz);

	parser reader_table(table);

	if(npart==5) {
//		ifstream inputfile("../src/SCEThaddec_3j/a1nnlofinite.txt");
		ifstream inputfile("../src/SCEThaddec_3j/a1nnlofinitelitreslog1.txt");
    	while (getline(inputfile, textread))
    	{
        	fline += textread;
    	}
	}
	else if(npart==10) {
//		ifstream inputfile("../src/SCEThaddec_3j/a2nnlofinite.txt");
		ifstream inputfile("../src/SCEThaddec_3j/a2nnlofinitelitreslog1.txt");
    	while (getline(inputfile, textread))
    	{
        	fline += textread;
    	}
	};

	ex expr = reader_table(fline).evalf();

	return ex_to<numeric>(expr).to_double();

}


