// -*- C++ -*-

#ifndef TPC_PAD_HELPER_HH
#define TPC_PAD_HELPER_HH

#include "DetectorID.hh"
static const int max_padid = 5768;
namespace tpc
{
const Double_t ZTarget = -143.; // Target from center
	enum EPadParameter
	{
		kLayerID,
		kNumOfPad,
		kRadius,
		kNumOfDivision,
		kDummy,
		kLength,
		NPadParameter
	};

	//_____________________________________________________________________________
	//#OfPad #division #radius padLength
	static const Double_t padParameter[NumOfLayersTPC][NPadParameter] =
	{{0, 48,    14.75, 48, 0,  9.},
		{1, 48,    24.25, 48, 0,  9.},
		{2, 72,    33.75, 72, 0,  9.},
		{3, 96,    43.25, 96, 0,  9.},
		{4, 120,    52.75,120,0,   9.},
		{5, 144,    62.25,144,0,   9.},
		{6, 168,    71.75,168,0,   9.},
		{7, 192,    81.25,192,0,   9.},
		{8, 216,    90.75,216,0,   9.},
		{9, 240,    100.25,240,0,  9.},
		{10,208,    111.5,241, 0,  12.5},
		{11,218,    124.5,271, 0,  12.5},
		{12,230,    137.5,300, 0,  12.5},
		{13,214,    150.5,330, 0,  12.5},
		{14,212,    163.5,360, 0,  12.5},
		{15,214,    176.5,390, 0,  12.5},
		{16,220,    189.5,420, 0,  12.5},
		{17,224,    202.5,449, 0,  12.5},
		{18,232,    215.5,479, 0,  12.5},
		{19,238,    228.5,509, 0,  12.5},
		{20,244,    241.5,539, 0,  12.5},
		{21,232,    254.5,569, 0,  12.5},
		{22,218,    267.5,599, 0,  12.5},
		{23,210,    280.5,628, 0,  12.5},
		{24,206,    293.5,658, 0,  12.5},
		{25,202,    306.5,688, 0,  12.5},
		{26,200,    319.5,718, 0,  12.5},
		{27,196,    332.5,748, 0,  12.5},
		{28,178,    345.5,777, 0,  12.5},
		{29,130,    358.5,807, 0,  12.5},
		{30,108,    371.5,837, 0,  12.5},
		{31,90,     384.5,867, 0, 12.5}};

	//_____________________________________________________________________________
	inline Int_t GetCoBoId(Int_t layer, Int_t row)
	{
		switch(layer){
			case 4:
				if(60<=row && row<=99)
					return 1;
				else
					return 0;
			case 5:
				if(72<=row && row<=119)
					return 1;
				else
					return 0;
			default:
				return layer/4;
		}
	}
static const Int_t padOnFrame[] =
{

//Pads on the frame
965,966,967,968,969,970,971,972,973,974,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1177,1178,1179,1180,1181,1182,1183,1184,1185,1186,1187,1188,1189,1190,1191,1192,1193,1194,1195,1196,1197,1198,1199,1200,1201,1202,1203,1204,1205,1206,1207,1208,1209,1210,1211,1212,1237,1238,1239,1240,1241,1242,1243,1244,1245,1246,1247,1248,1249,1250,1251,1252,1253,1254,1255,1256,1257,1258,1259,1260,1261,1262,1263,1264,1265,1266,1267,1268,1269,1270,1271,1272,1395,1396,1397,1398,1399,1400,1401,1405,1406,1407,1408,1409,1410,1411,1412,1413,1414,1415,1416,1417,1418,1419,1420,1421,1422,1423,1424,1425,1426,1427,1428,1429,1430,1431,1432,1436,1437,1438,1439,1440,1441,1442,1455,1456,1457,1458,1459,1460,1461,1465,1466,1467,1468,1469,1470,1471,1472,1473,1474,1475,1476,1477,1478,1479,1480,1481,1482,1483,1484,1485,1486,1487,1488,1489,1490,1491,1492,1496,1497,1498,1499,1500,1501,1502,1595,1596,1597,1598,1604,1605,1606,1607,1608,1647,1648,1649,1650,1651,1652,1658,1659,1660,1663,1664,1665,1671,1672,1673,1674,1675,1676,1715,1716,1717,1718,1719,1725,1726,1727,1728,1807,1808,1809,1810,1811,1812,1813,1814,1815,1816,1817,1818,1879,1880,1881,1882,1889,1890,1891,1892,1953,1954,1955,1956,1957,1958,1959,1960,1961,1962,1963,1964,2018,2019,2020,2025,2026,2027,2028,2100,2101,2102,2103,2104,2111,2112,2113,2114,2115,2187,2188,2189,2190,2195,2196,2197,2219,2220,2221,2222,2223,2224,2225,2226,2227,2228,2309,2310,2311,2317,2318,2323,2324,2330,2331,2332,2413,2414,2415,2416,2417,2418,2419,2420,2421,2422,2427,2428,2429,2430,2518,2519,2520,2521,2522,2523,2524,2525,2526,2527,2540,2541,2542,2543,2544,2545,2546,2547,2548,2549,2637,2638,2639,2640,2732,2733,2739,2740,2761,2762,2768,2769,2950,2951,2952,2953,2954,2955,2956,2957,2958,2987,2988,2989,2990,2991,2992,2993,2994,2995,3174,3175,3176,3177,3178,3179,3180,3181,3182,3219,3220,3221,3222,3223,3224,3225,3226,3227,3405,3406,3407,3408,3409,3410,3411,3412,3413,3458,3459,3460,3461,3462,3463,3464,3465,3466,3642,3643,3644,3645,3646,3647,3648,3649,3650,3703,3704,3705,3706,3707,3708,3709,3710,3711,3877,3878,3879,3880,3881,3882,3883,3884,3945,3946,3947,3948,3949,3950,3951,3952,4098,4099,4105,4174,4180,4181,4308,4309,4310,4311,4312,4313,4314,4315,4316,4391,4392,4393,4394,4395,4396,4397,4398,4399,4512,4513,4520,4603,4610,4611,4713,4714,4715,4716,4717,4718,4719,4720,4811,4812,4813,4814,4815,4816,4817,4818,4910,4911,4912,4913,4914,4915,4916,4917,5016,5017,5018,5019,5020,5021,5022,5023,5104,5105,5108,5111,5112,5217,5218,5221,5224,5225,5287,5288,5289,5290,5291,5292,5293,5294,5295,5408,5409,5410,5411,5412,5413,5414,5415,5416,5441,5442,5443,5444,5445,5566,5567,5568,5569,5570,

//Empty Pads
1394,1402,1403,1404,1433,1434,1435,1462,1463,1464,1493,1494,1495,1503,1599,1601,1603,1653,1655,1657,1666,1667,1668,1669,1670,1720,1721,1722,1723,1724,1877,1878,1883,1884,1885,1886,1887,1888,1893,1894,2021,2022,2023,2024,2105,2106,2107,2108,2109,2110,2191,2192,2193,2194,2312,2313,2314,2315,2316,2325,2326,2327,2328,2329,2734,2735,2736,2737,2738,2748,2763,2764,2765,2766,2767,4100,4101,4102,4103,4104,4175,4176,4177,4178,4179,4514,4515,4516,4517,4518,4519,4604,4605,4606,4607,4608,4609,5106,5107,5109,5110,5219,5220,5222,5223
};

	//_____________________________________________________________________________
	inline Int_t GetPadId(Int_t layerID, Int_t rowID)
	{
		//Original
		//Int_t padID=0;
		// Check!!!!!!!
		Int_t padID=1;
		for(int layi = 0 ; layi<layerID; layi++) padID += padParameter[layi][1];
		padID+=rowID;
		return padID;

	}

	inline Int_t getLayerID(Int_t padID)
	{
		padID-=1;
		int layer;
		int sum = 0;

		for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
		{
			sum += padParameter[layer][1];
		}
		return layer;
	}

	inline Int_t getRowID(Int_t padID)
	{
		padID-=1;
		int layer, row;
		int sum = 0;

		for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
		{
			sum += padParameter[layer][1];
		}
		row = padID - sum;
		return row;
	}
	/*
		 Double_t getTheta(Int_t layerID, Int_t rowID)
		 {
		 Double_t sTheta = 180.-(360./padParameter[layerID][3])*padParameter[layerID][1]/2.;
		 Double_t theta = sTheta+(rowID+0.5)*(360.-2*sTheta)/padParameter[layerID][1];
		 return theta;
		 }
		 */

	inline Double_t getTheta(Int_t padID)
	{
		padID-=1;
		int layer, row;
		int sum = 0;

		for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
		{
			sum += padParameter[layer][1];
		}
		row = padID - sum;
		//std::cout<<"layer="<<layer<<", row="<<row<<std::endl;
		Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;
		//Double_t theta = sTheta+(row+0.5)*(360.-2*sTheta)/padParameter[layer][1];
		Double_t theta = sTheta+(row+0.5)*360./padParameter[layer][3]-180;
		//std::cout<<"theta="<<theta<<std::endl;

		return theta;
	}

	inline Double_t getTheta(Int_t layer, Double_t m_row)
	{
		Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;
		//Double_t theta = sTheta+(row+0.5)*(360.-2*sTheta)/padParameter[layer][1];
		Double_t theta = sTheta+(m_row+0.5)*360./padParameter[layer][3]-180;
		//std::cout<<"theta="<<theta<<std::endl;
		return theta;
	}


	inline Double_t getR(Int_t padID)
	{
		padID-=1;
		int layer;
		int sum = 0;

		for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
		{
			sum += padParameter[layer][1];
		}
		Double_t R = padParameter[layer][2];
		return R;
	}

	inline TVector3 getPosition(int padID)
	{
		padID-=1;
		int layer, row;
		int sum = 0;

		for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
		{
			sum += padParameter[layer][1];
		}
		row = padID - sum;

		TVector3 result;
		if (row > padParameter[layer][1]){ // out of range
			result.SetX(0);
			result.SetY(-1);
			result.SetZ(0);
		}
		else{
			double x, z;
			//Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;

			//    x = padParameter[layer][2] * -sin((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.);
			//    z = padParameter[layer][2] * -cos((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.) - 143.0;

			// x = padParameter[layer][2] * sin(getTheta(padID+1)*TMath::Pi()/180.);
			// z = padParameter[layer][2] * cos(getTheta(padID+1)*TMath::Pi()/180.) - 143.0;
			//std::cout<<"layer="<<layer<<", row"<<row<<std::endl;
			x = padParameter[layer][2] * sin(getTheta(layer,row)*TMath::Pi()/180.);
			z = padParameter[layer][2] * cos(getTheta(layer,row)*TMath::Pi()/180.) - 143.0;

			// Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;
			// double x_ = padParameter[layer][2] * -sin((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.);
			// double z_ = padParameter[layer][2] * -cos((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.) - 143.0;
			//std::cout<<"x="<<x<<", z="<<z<<", x_="<<x_<<", z_="<<z_<<std::endl;
			result.SetX(x);
			result.SetY(0);
			result.SetZ(z);
		}
		return result;
	}

	inline TVector3 getPosition(int layer, double m_row)
	{
		TVector3 result;
		if (m_row > padParameter[layer][1]){ // out of range
			result.SetX(0);
			result.SetY(-1);
			result.SetZ(0);
		}
		else{
			double x, z;

			x = padParameter[layer][2] * sin(getTheta(layer, m_row)*TMath::Pi()/180.);
			z = padParameter[layer][2] * cos(getTheta(layer, m_row)*TMath::Pi()/180.) - 143.0;
			result.SetX(x);
			result.SetY(0);
			result.SetZ(z);
			//std::cout<<"x="<<x<<", z="<<z<<std::endl;
		}
		return result;
	}



	inline int findPadID(double z, double x)
	{
		z += 143;
		double radius = sqrt(x*x + z*z);
		double angle;
		if (z == 0)
		{
			if (x > 0)   angle = 1.5*TMath::Pi();
			else if (x < 0)   angle = 0.5*TMath::Pi();
			else return -1000; // no padID if (0,0)
		}
		//  else
		//  {
		//    if (z < 0 && x < 0) angle = atan(x / z);
		//    else if (z > 0 && x < 0) angle = TMath::Pi() - atan(-x / z);
		//    else if (z > 0 && x > 0) angle = TMath::Pi() + atan(x / z);
		//    else if (z < 0 && x > 0) angle = 2*TMath::Pi() - atan(-x / z);
		//  }
		else if (z > 0) angle = TMath::Pi()+atan(x / z);
		else if( z < 0&&x<0) angle = atan(x / z);
		else angle = 2*TMath::Pi()+ atan(x / z);
		//cout << " angle: " << angle*180/TMath::Pi() << endl;
		if(z<0&&x>0){
			//		cout<<Form("Angle: %f,  (%f , %f ) ",angle,z,x)<<endl;
		}
		int layer, row;
		// find layer_num.
		for (layer = 0; !(padParameter[layer][2]+padParameter[layer][5]*0.5 >= radius
					&& padParameter[layer][2]-padParameter[layer][5]*0.5 <= radius); layer++)
		{
			if (layer >= 32) return -1000;
			if (layer != 0)
			{
				if (padParameter[layer][2] - padParameter[layer][5] * 0.5 >= radius &&
						padParameter[layer - 1][2] + padParameter[layer - 1][5] * 0.5 <= radius) return -layer;
			}
		}

		//cout << " layer: " << layer << endl;

		Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;

		// find row_num
		if (angle - (sTheta*TMath::Pi()/180.) < 0)
		{
//			cout<<Form("Angle: %f,  (%f , %f ) ",angle,z,x)<<endl;
			return -1000;
		}

		//double a, b, c;
		row = (int)((angle-(sTheta*TMath::Pi()/180.))/(360./padParameter[layer][3]*TMath::Pi()/180.));
		if (row > padParameter[layer][1]) return -1000;

		//cout << " row: " << row << endl;
		//This is original one
		//return GetPadId(layer, row)+1;
		//Please check
		return GetPadId(layer, row);
	}

	//_____________________________________________________________________________
	inline
		TH2Poly* InitializeHistograms()
		{
			/*
				 std::vector<TH2Poly*> target;
				 TList* list = gDirectory->GetList();
				 TIter itr(list);
				 while(itr.Next()){
				 const TString& name((*itr)->GetName());
				 const TString& cname((*itr)->ClassName());
			//hddaq::cout << " " << std::setw(8) << std::left << name
			std::cout << " " << std::setw(8) << std::left << name
			<< "(" << cname << ")" << std::endl;
			if(cname.EqualTo("TH2Poly")){
			target.push_back(dynamic_cast<TH2Poly*>(*itr));
			}
			}
			*/
			Double_t X[5];
			Double_t Y[5];
			TH2Poly* h = new TH2Poly();
			for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
				Double_t pLength = tpc::padParameter[layer][5];
				Double_t st = 180.-(360./tpc::padParameter[layer][3])
					* tpc::padParameter[layer][1]/2.;
				Double_t sTheta  = (-1+st/180.)*TMath::Pi();
				Double_t dTheta  = (360./tpc::padParameter[layer][3])/180.*TMath::Pi();
				Double_t cRad    = tpc::padParameter[layer][2];
				Int_t    nPad    = tpc::padParameter[layer][1];
				for(Int_t j=0; j<nPad; ++j){
					X[1] = (cRad+(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
					X[2] = (cRad+(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
					X[3] = (cRad-(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
					X[4] = (cRad-(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
					X[0] = X[4];
					Y[1] = (cRad+(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
					Y[2] = (cRad+(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
					Y[3] = (cRad-(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
					Y[4] = (cRad-(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
					Y[0] = Y[4];
					for(Int_t ii=0; ii<5; ++ii) X[ii] -=143;
					//      for(auto& h: target) h->AddBin(5, X, Y);
					h->AddBin(5, X, Y);
				}
			}
			return h;
		}
}
void GetLayerRow(int padId, int &layer, int &row){
	layer=tpc::getLayerID(padId);
	row=tpc::getRowID(padId);
};
#endif
