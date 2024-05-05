#pragma once
#define MAXSATNUM 100
#define MAXM (MAXSATNUM*4)
#define MAXN (MAXSATNUM*2+3)

#define c_light 2.99792458E8//光速
#define FG1_GPS  1575.42E6             /* GPS L1信号频率 */
#define FG2_GPS  1227.60E6             /* GPS L2信号频率 */
#define WL1_GPS  (c_light/FG1_GPS)
#define WL2_GPS  (c_light/FG2_GPS)

#define BDS_GM 3.986004418E14//BDS系统的GM
#define BDS_OMEGA 7.2921150E-5//BDS系统地球自转角速度
#define FG1_BDS  1561.098E6               /* B1信号的基准频率 */
#define FG2_BDS  1207.140E6               /* B2信号的基准频率 */
#define FG3_BDS  1268.520E6               /* B2信号的基准频率 */
#define WL1_BDS  (c_light/FG1_BDS)
#define WL2_BDS  (c_light/FG2_BDS)
#define WL3_BDS  (c_light/FG3_BDS)
enum System {GPS,BDS,UNKNOWN};
struct Time {
	int week;
	double sow;
};
struct SatPos {
	System sys;
	int prn;
	bool is_refsat;
	double posrov[3];
	double posbas[3];
};
struct Obs{
	System sys;
	int prn;
	double P[2];
	double L[2];
};
struct Data {
	Time time;
	int satnum[2];
	int refsatidx[2];
	Obs obs[MAXSATNUM];
	SatPos satpos[MAXSATNUM];
};
