#ifndef __CONSTANTES_H__
#define __CONSTANTES_H__
//Datos generales
const int num_tipoPlato = 3;
const int num_nutr = 21;
const int num_alerg = 7;
const int num_incomp = 5;
const int num_gAl = 10;
const int num_obj = 2;

//Cantidad de nutrientes recomedados (por almuerzo)
const int CALCIUM_INDEX = 1;
const int ENERGY_INDEX = 2;
const int FAT_INDEX = 4;
const int IRON_INDEX = 5;
const int POTASIUM_INDEX = 7;
const int PROTEIN_INDEX = 8;


const int FORCED_INDEXES_SIZE = 3;
const int FORCED_INDEXES[3] = {ENERGY_INDEX, FAT_INDEX, PROTEIN_INDEX};
const double FORCED_MIN[3] = {0.85, 0.75, 0.75};
const double FORCED_MAX[3] = {1.15, 1.25, 1.25};

const double ingR_acFol = 135;  //acido folico
const double ingR_cal = 585;    //calcio
const double ingR_en = 1012;    //energia
const double ingR_fos = 562.5;  //fosforo
const double ingR_gra = 31.72;  //grasa
const double ingR_hie = 8.55;   //hierro
const double ingR_mag = 112.5;  //magnesio
const double ingR_pot = 2025;   //potasio
const double ingR_pro = 27.08;  //proteinas
const double ingR_sel = 25.75;  //selenio
const double ingR_sod = 870;    //sodio
const double ingR_vA = 450;     //vitA
const double ingR_vB1 = 0.41;   //vitB1
const double ingR_vB2 = 0.63;   //vitB2
const double ingR_vB6 = 0.54;   //vitB6
const double ingR_vB12 = 2.28;  //vitB12
const double ingR_vC = 27;      //vitC
const double ingR_vD = 4.65;    //vitD
const double ingR_vE = 6.3;     //vitE
const double ingR_yod = 67.5;   //yodo
const double ingR_zin = 6.75;   //zinc

const double ingR[21] = {ingR_acFol, ingR_cal, ingR_en, ingR_fos, ingR_gra, ingR_hie, ingR_mag, ingR_pot, ingR_pro, ingR_sel, ingR_sod, ingR_vA, ingR_vB1, ingR_vB2, ingR_vB6, ingR_vB12, ingR_vC, ingR_vD, ingR_vE, ingR_yod, ingR_zin};
const double minReq[21] = {0.70, 0.70, 0.90, 0.70, 0.90, 0.70, 0.70, 0.70, 0.90, 0.70, 0.70, 0.70, 0.70, 0.70, 0.70, 0.70, 0.70, 0.70, 0.70, 0.70, 0.70};
const double maxReq[21] = {1.30, 1.30, 1.10, 1.30, 1.10, 1.30, 1.30, 1.30, 1.10, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.30};

const int NPLATOS[3] = {19, 34, 14};
#endif
