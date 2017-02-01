#include "TimeIntegration.h"

#include <iostream>

using namespace std;

void criaParticulasQuebraBarragem2D(){

	int mp = 36;
	int np = 36;

	double xl = mp*0.008;
	double yl = np*0.008;

	double dx = xl / mp;
	double dy = yl / np;

	double vx = 0.0;

	FILE *writeGrid;
	writeGrid = fopen("mps.grd", "w");
	//fprintf(writeGrid, "0.0\n");
	fprintf(writeGrid, "%d \n", ((mp / 2)*np) + (4 * mp + 3 * 1 + 7) + (4 * mp + 3 * 3 + 5) + (4 * mp + 3 * 5 + 3) /*+ (4 * mp + 7 + 8 + 7) + 3 * 8*/);

	//PARTICULAS TIPO 0 #fluido
	for (int i = 0; i < (mp / 2); i++){
		for (int j = 0; j < np; j++){
			fprintf(writeGrid, "0 %lf %lf 0.0 0.0 0.0 0.0 0.0 0.0 \n", (i*dx + (3 * dx)), (j*dy + (3 * dy)));
		}
	}

	//PARTICULAS TIPO 2 #interna
	//for (int i = 0; i < (mp) + 1; ++i){ //up
	//	fprintf(writeGrid, "2  %lf %lf %lf 0.0 0.0 0.0 \n", ((i*dx) + 3*dx), (yl + 3*dy), vx);
	//}
	for (int i = 0; i < (2 * mp) + 1; ++i){ //down
		fprintf(writeGrid, "2 %lf %lf 0.0 0.0 0.0 0.0 0.0 0.0 \n", ((i*dx) + 2 * dx), (2 * dy));
	}
	for (int i = 0; i < (mp)+4; ++i){ //left
		fprintf(writeGrid, "2 %lf %lf 0.0 0.0 0.0 0.0 0.0 0.0 \n", (2 * dx), ((i*dy) + (3 * dy)));
	}
	for (int i = 0; i < (mp)+5; ++i){ //right
		fprintf(writeGrid, "2 %lf %lf 0.0 0.0 0.0 0.0 0.0 0.0 \n", (2 * xl + 3 * dx), ((i*dy) + (2 * dy)));
	}

	//PARTICULAS TIPO 2 #externa
	//for (int i = 0; i < (mp) + 3; ++i){ //up
	//  fprintf(writeGrid, "3  %lf %lf %lf 0.0 0.0 0.0 \n", ((i*dx) + (2*dx)), (yl + (4*dy)), vx);
	//}
	for (int i = 0; i < (2 * mp) + 3; ++i){ //down
		fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", ((i*dx) + dx), dy, vx);
	}
	for (int i = 0; i < (mp)+5; ++i){ //left
		fprintf(writeGrid, "3 %lf %lf 0.0 0.0 0.0 0.0 0.0 0.0 \n", dx, ((i*dy) + (2 * dy)));
	}
	for (int i = 0; i < (mp)+6; ++i){ //right
		fprintf(writeGrid, "3 %lf %lf 0.0 0.0 0.0 0.0 0.0 0.0 \n", (2 * xl + (4 * dx)), ((i*dy) + dy));
	}

	//PARTICULAS TIPO 3 
	//for (int i = 0; i < mp + 5; ++i){ //up
	//	fprintf(writeGrid,  ""3 " << i*dx + dx << " " << yl + (5*dy) << " " << vx << " 0.0 0.0 0.0" << endl;
	//}
	for (int i = 0; i < (2 * mp) + 5; ++i){ //down
		fprintf(writeGrid, "3 %lf 0.0 %lf 0.0 0.0 0.0 0.0 0.0 \n", (i*dx), vx);
	}
	for (int i = 0; i < (mp)+6; ++i){ //left
		fprintf(writeGrid, "3 0.0 %lf 0.0 0.0 0.0 0.0 0.0 0.0 \n", ((i*dy) + dy));
	}
	for (int i = 0; i < (mp)+7; ++i){ //right
		fprintf(writeGrid, "3 %lf %lf 0.0 0.0 0.0 0.0 0.0 0.0 \n", (2 * xl + (5 * dx)), (i*dy));
	}

	//PARTICULAS TIPO 3 #externa
	////for (int i = 0; i < mp + 5; ++i){ //up
	////	fprintf(writeGrid,  ""3 " << i*dx + dx << " " << yl + (5*dy) << " " << vx << " 0.0 0.0 0.0" << endl;
	////}
	//for (int i = 0; i < (2 * mp) + 7; ++i){ //down
	//	fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 \n", (i*dx), -dy, vx);
	//}
	//for (int i = 0; i < (mp)+8; ++i){ //left
	//	fprintf(writeGrid, "3 %lf %lf 0.0 0.0 0.0 0.0 \n", -dx, (i*dy) - dy);
	//}
	//for (int i = 0; i < (mp)+7; ++i){ //right
	//	fprintf(writeGrid, "3 %lf %lf 0.0 0.0 0.0 0.0 \n", (2 * xl + (6 * dx)), (i*dy));
	//}

	//PARTICULAS DO TOCO DO TIPO 2
	/*for (int i = 0; i < 3; ++i){
	for (int j = 0; j < 8; j++){
	fprintf(writeGrid, "2 %lf %lf 0.0 0.0 0.0 0.0 \n", ((i*dx) + (mp + 1) * dx), (3 * dy) + j*dy);
	}
	}*/

	fclose(writeGrid);
}

void criaParticulasQuebraBarragem3D(){

	int mp = 14;
	int np = 7;
	int op = 7;

	double xl = mp*0.0125;
	double yl = np*0.0125;
	double zl = op*0.0125;

	double dx = xl / mp;
	double dy = yl / np;
	double dz = zl / op;

	//double vx = 0.0;

	FILE *writeGrid;
	writeGrid = fopen("mps.grd", "w");
	//fprintf(writeGrid, "0.0\n");
	fprintf(writeGrid, "%d \n", /*fluido*/(op)*(mp*np) +/*first wall*/op*(4 * mp + 10) + (4 * mp + 4)*(mp + 5) +/*second wall*/(op + 4)*(4 * mp + 14) + (4 * mp + 4)*(mp + 5) +/*thrid wall*/(op + 6)*(4 * mp + 18) + (4 * mp + 8)*(mp + 6) /*+ /*fourth wall (op + 8)*(4 * mp + 22) + (4 * mp + 12)*(mp + 7) + /*fifth wall (op + 10)*(4 * mp + 26) + (4 * mp + 16)*(mp + 8)*/);

	for (int k = 0; k < (op); k++){
		//PARTICULAS TIPO 0 #fluido
		for (int i = 0; i < (mp / 2); i++){
			for (int j = 0; j < (2 * np); j++){
				fprintf(writeGrid, "0 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", (i*dx + (3 * dx)), (j*dy + (3 * dy)), dz*k);
			}
		}
	}

	//PARTICULAS TIPO 2 #interna
	for (int k = 0; k < (op); k++){
		for (int i = 0; i < (2 * mp) + 1; ++i){ //down
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", ((i*dx) + 2 * dx), (2 * dy), dz*k);
		}
		for (int i = 0; i < (mp)+4; ++i){ //left
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", (2 * dx), ((i*dy) + (3 * dy)), dz*k);
		}
		for (int i = 0; i < (mp)+5; ++i){ //right
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", (2 * xl + 3 * dx), ((i*dy) + (2 * dy)), dz*k);
		}
	}
	for (int i = 0; i < ((2 * mp) + 2); ++i){ //front
		for (int j = 0; j < ((mp)+5); ++j){
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", ((i*dx) + 2 * dx), ((2 * dy) + j*dy), -dz);
		}
	}
	for (int i = 0; i < ((2 * mp) + 2); ++i){ //back
		for (int j = 0; j < ((mp)+5); ++j){
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", ((i*dx) + 2 * dx), ((2 * dy) + j*dy), (op)*dz);
		}
	}

	//PARTICULAS TIPO 2 #externa
	for (int k = 0; k < (op + 4); k++){
		for (int i = 0; i < (2 * mp) + 4; ++i){ //down
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", ((i*dx) + dx), dy, dz*k - 2 * dz);
		}
		for (int i = 0; i < (mp)+5; ++i){ //left
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", dx, ((i*dy) + (2 * dy)), dz*k - 2 * dz);
		}
		for (int i = 0; i < (mp)+5; ++i){ //right
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", (2 * xl + (4 * dx)), ((i*dy) + 2 * dy), dz*k - 2 * dz);
		}
	}

	for (int i = 0; i < ((2 * mp) + 2); ++i){ //front
		for (int j = 0; j < ((mp)+5); ++j){
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", ((i*dx) + 2 * dx), (2 * dy + j*dy), -2 * dz);
		}
	}
	for (int i = 0; i < ((2 * mp) + 2); ++i){ //back
		for (int j = 0; j < ((mp)+5); ++j){
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", ((i*dx) + 2 * dx), (2 * dy + j*dy), (op + 1)*dz);
		}
	}

	//PARTICULAS TIPO 3 #interna
	for (int k = 0; k < (op + 6); k++){
		for (int i = 0; i < (2 * mp) + 6; ++i){ //down
			fprintf(writeGrid, "3 %lf 0.0 %lf 0.0 0.0 0.0 0.0 0.0 \n", (i*dx), dz*k - 3 * dz);
		}
		for (int i = 0; i < (mp)+6; ++i){ //left
			fprintf(writeGrid, "3 0.0 %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", ((i*dy) + dy), dz*k - 3 * dz);
		}
		for (int i = 0; i < (mp)+6; ++i){ //right
			fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", (2 * xl + (5 * dx)), (i*dy + dy), dz*k - 3 * dz);
		}
	}

	for (int i = 0; i < ((2 * mp) + 4); ++i){ //front
		for (int j = 0; j < ((mp)+6); ++j){
			fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", i*dx + dx, j*dy + dy, -3 * dz);
		}
	}
	for (int i = 0; i < ((2 * mp) + 4); ++i){ //back
		for (int j = 0; j < ((mp)+6); ++j){
			fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", i*dx + dx, j*dy + dy, (op + 2)*dz);
		}
	}

	///PARTICULAS TIPO 3 #meio
	/*for (int k = 0; k < (op + 8); k++){
		for (int i = 0; i < (2 * mp) + 8; ++i){ //down
			fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", ((i*dx) - dx), (-1 * dy), (-4 * dz) + dz*k);
		}
		for (int i = 0; i < (mp)+7; ++i){ //left
			fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", (-dx), (i*dy), (-4 * dz) + dz*k);
		}
		for (int i = 0; i < (mp)+7; ++i){ //right
			fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", (2 * xl + 6 * dx), (i*dy), (-4 * dz) + dz*k);
		}
	}
	for (int i = 0; i < ((2 * mp) + 6); ++i){ //front
		for (int j = 0; j < ((mp)+7); ++j){
			fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", (i*dx), (j*dy), -4 * dz);
		}
	}
	for (int i = 0; i < ((2 * mp) + 6); ++i){ //back
		for (int j = 0; j < ((mp)+7); ++j){
			fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", (i*dx), (j*dy), (op + 3)*dz);
		}
	}*/

	///PARTICULAS TIPO 3 #externa
	/*for (int k = 0; k < (op + 10); k++){
		for (int i = 0; i < (2 * mp) + 10; ++i){ //down
			fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", ((i*dx) - 2*dx), (-2 * dy), (-5 * dz) + dz*k);
		}
		for (int i = 0; i < (mp)+8; ++i){ //left
			fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", (-2*dx), (i*dy)-dy, (-5 * dz) + dz*k);
		}
		for (int i = 0; i < (mp)+8; ++i){ //right
			fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", (2 * xl + 7 * dx), (i*dy)-dy, (-5 * dz) + dz*k);
		}
	}
	for (int i = 0; i < ((2 * mp) + 8); ++i){ //front
		for (int j = 0; j < ((mp)+8); ++j){
			fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", (i*dx)-dx, (j*dy)-dy, -5 * dz);
		}
	}
	for (int i = 0; i < ((2 * mp) + 8); ++i){ //back
		for (int j = 0; j < ((mp)+8); ++j){
			fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", (i*dx)-dx, (j*dy)-dy, (op + 4)*dz);
		}
	}*/

	fclose(writeGrid);
}

void criaParticulasQuebraBarragem3D_offset() {

	int mp = 14;
	int np = 7;
	int op = 7;

	double xl = mp*0.0125;
	double yl = np*0.0125;
	double zl = op*0.0125;

	double dx = xl / mp;
	double dy = yl / np;
	double dz = zl / op;

	double offset = 3*dx;

	FILE *writeGrid;
	writeGrid = fopen("mps.grd", "w");
	//fprintf(writeGrid, "0.0\n");
	fprintf(writeGrid, "%d \n", /*fluido*/(op)*(mp*np) +/*first wall*/op*(4 * mp + 10) + (4 * mp + 4)*(mp + 5) +/*second wall*/(op + 4)*(4 * mp + 14) + (4 * mp + 4)*(mp + 5) +/*thrid wall*/(op + 6)*(4 * mp + 18) + (4 * mp + 8)*(mp + 6) /*+ /*fourth wall (op + 8)*(4 * mp + 22) + (4 * mp + 12)*(mp + 7) + /*fifth wall (op + 10)*(4 * mp + 26) + (4 * mp + 16)*(mp + 8)*/);

	for (int k = 0; k < (op); k++) {
		//PARTICULAS TIPO 0 #fluido
		for (int i = 0; i < (mp / 2); i++) {
			for (int j = 0; j < (2 * np); j++) {
				fprintf(writeGrid, "0 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", (i*dx + (3 * dx)) + offset, (j*dy + (3 * dy)) + offset, dz*k + offset);
			}
		}
	}

	//PARTICULAS TIPO 2 #interna
	for (int k = 0; k < (op); k++) {
		for (int i = 0; i < (2 * mp) + 1; ++i) { //down
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", ((i*dx) + 2 * dx) + offset, (2 * dy) + offset, dz*k + offset);
		}
		for (int i = 0; i < (mp)+4; ++i) { //left
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", (2 * dx) + offset, ((i*dy) + (3 * dy)) + offset, dz*k + offset);
		}
		for (int i = 0; i < (mp)+5; ++i) { //right
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", (2 * xl + 3 * dx) + offset, ((i*dy) + (2 * dy)) + offset, dz*k + offset);
		}
	}
	for (int i = 0; i < ((2 * mp) + 2); ++i) { //front
		for (int j = 0; j < ((mp)+5); ++j) {
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", ((i*dx) + 2 * dx) + offset, ((2 * dy) + j*dy) + offset, -dz + offset);
		}
	}
	for (int i = 0; i < ((2 * mp) + 2); ++i) { //back
		for (int j = 0; j < ((mp)+5); ++j) {
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", ((i*dx) + 2 * dx) + offset, ((2 * dy) + j*dy) + offset, (op)*dz + offset);
		}
	}

	//PARTICULAS TIPO 2 #externa
	for (int k = 0; k < (op + 4); k++) {
		for (int i = 0; i < (2 * mp) + 4; ++i) { //down
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", ((i*dx) + dx) + offset, dy + offset, dz*k - 2 * dz + offset);
		}
		for (int i = 0; i < (mp)+5; ++i) { //left
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", dx + offset, ((i*dy) + (2 * dy)) + offset, dz*k - 2 * dz + offset);
		}
		for (int i = 0; i < (mp)+5; ++i) { //right
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", (2 * xl + (4 * dx)) + offset, ((i*dy) + 2 * dy) + offset, dz*k - 2 * dz + offset);
		}
	}

	for (int i = 0; i < ((2 * mp) + 2); ++i) { //front
		for (int j = 0; j < ((mp)+5); ++j) {
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", ((i*dx) + 2 * dx) + offset, (2 * dy + j*dy) + offset, -2 * dz + offset);
		}
	}
	for (int i = 0; i < ((2 * mp) + 2); ++i) { //back
		for (int j = 0; j < ((mp)+5); ++j) {
			fprintf(writeGrid, "2 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", ((i*dx) + 2 * dx) + offset, (2 * dy + j*dy) + offset, (op + 1)*dz + offset);
		}
	}

	//PARTICULAS TIPO 3 #interna
	for (int k = 0; k < (op + 6); k++) {
		for (int i = 0; i < (2 * mp) + 6; ++i) { //down
			fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", (i*dx) + offset, offset,dz*k - 3 * dz + offset);
		}
		for (int i = 0; i < (mp)+6; ++i) { //left
			fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", offset, ((i*dy) + dy) + offset, dz*k - 3 * dz + offset);
		}
		for (int i = 0; i < (mp)+6; ++i) { //right
			fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", (2 * xl + (5 * dx)) + offset, (i*dy + dy) + offset, dz*k - 3 * dz + offset);
		}
	}

	for (int i = 0; i < ((2 * mp) + 4); ++i) { //front
		for (int j = 0; j < ((mp)+6); ++j) {
			fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", i*dx + dx + offset, j*dy + dy + offset, -3 * dz + offset);
		}
	}
	for (int i = 0; i < ((2 * mp) + 4); ++i) { //back
		for (int j = 0; j < ((mp)+6); ++j) {
			fprintf(writeGrid, "3 %lf %lf %lf 0.0 0.0 0.0 0.0 0.0 \n", i*dx + dx + offset, j*dy + dy + offset, (op + 2)*dz + offset);
		}
	}

	
	fclose(writeGrid);
}

void criaParticulasQuebraBarragem3D_offset_semzero() {

	int mp = 14;
	int np = 7;
	int op = 7;

	double xl = mp*0.0125;
	double yl = np*0.0125;
	double zl = op*0.0125;

	double dx = xl / mp;
	double dy = yl / np;
	double dz = zl / op;

	double offsetX = 2 * dx;
	double offsetY = 2 * dx;
	double offsetZ = 4 * dx;

	FILE *writeGrid;
	writeGrid = fopen("mps.grd", "w");
	//fprintf(writeGrid, "0.0\n");
	fprintf(writeGrid, "%d \n", /*fluido*/(op)*(mp*np) +/*first wall*/op*(4 * mp + 10) + (4 * mp + 4)*(mp + 5) +/*second wall*/(op + 4)*(4 * mp + 14) + (4 * mp + 4)*(mp + 5) +/*thrid wall*/(op + 6)*(4 * mp + 18) + (4 * mp + 8)*(mp + 6) /*+ /*fourth wall (op + 8)*(4 * mp + 22) + (4 * mp + 12)*(mp + 7) + /*fifth wall (op + 10)*(4 * mp + 26) + (4 * mp + 16)*(mp + 8)*/);

	for (int k = 0; k < (op); k++) {
		//PARTICULAS TIPO 0 #fluido
		for (int i = 0; i < (mp / 2); i++) {
			for (int j = 0; j < (2 * np); j++) {
				fprintf(writeGrid, "0 %lf %lf %lf\n", (i*dx + (3 * dx)) + offsetX, (j*dy + (3 * dy)) + offsetY, dz*k + offsetZ);
			}
		}
	}

	//PARTICULAS TIPO 2 #interna
	for (int k = 0; k < (op); k++) {
		for (int i = 0; i < (2 * mp) + 1; ++i) { //down
			fprintf(writeGrid, "2 %lf %lf %lf\n", ((i*dx) + 2 * dx) + offsetX, (2 * dy) + offsetY, dz*k + offsetZ);
		}
		for (int i = 0; i < (mp)+4; ++i) { //left
			fprintf(writeGrid, "2 %lf %lf %lf\n", (2 * dx) + offsetX, ((i*dy) + (3 * dy)) + offsetY, dz*k + offsetZ);
		}
		for (int i = 0; i < (mp)+5; ++i) { //right
			fprintf(writeGrid, "2 %lf %lf %lf\n", (2 * xl + 3 * dx) + offsetX, ((i*dy) + (2 * dy)) + offsetY, dz*k + offsetZ);
		}
	}
	for (int i = 0; i < ((2 * mp) + 2); ++i) { //front
		for (int j = 0; j < ((mp)+5); ++j) {
			fprintf(writeGrid, "2 %lf %lf %lf\n", ((i*dx) + 2 * dx) + offsetX, ((2 * dy) + j*dy) + offsetY, -dz + offsetZ);
		}
	}
	for (int i = 0; i < ((2 * mp) + 2); ++i) { //back
		for (int j = 0; j < ((mp)+5); ++j) {
			fprintf(writeGrid, "2 %lf %lf %lf\n", ((i*dx) + 2 * dx) + offsetX, ((2 * dy) + j*dy) + offsetY, (op)*dz + offsetZ);
		}
	}

	//PARTICULAS TIPO 2 #externa
	for (int k = 0; k < (op + 4); k++) {
		for (int i = 0; i < (2 * mp) + 4; ++i) { //down
			fprintf(writeGrid, "2 %lf %lf %lf\n", ((i*dx) + dx) + offsetX, dy + offsetY, dz*k - 2 * dz + offsetZ);
		}
		for (int i = 0; i < (mp)+5; ++i) { //left
			fprintf(writeGrid, "2 %lf %lf %lf\n", dx + offsetX, ((i*dy) + (2 * dy)) + offsetY, dz*k - 2 * dz + offsetZ);
		}
		for (int i = 0; i < (mp)+5; ++i) { //right
			fprintf(writeGrid, "2 %lf %lf %lf\n", (2 * xl + (4 * dx)) + offsetX, ((i*dy) + 2 * dy) + offsetY, dz*k - 2 * dz + offsetZ);
		}
	}

	for (int i = 0; i < ((2 * mp) + 2); ++i) { //front
		for (int j = 0; j < ((mp)+5); ++j) {
			fprintf(writeGrid, "2 %lf %lf %lf\n", ((i*dx) + 2 * dx) + offsetX, (2 * dy + j*dy) + offsetY, -2 * dz + offsetZ);
		}
	}
	for (int i = 0; i < ((2 * mp) + 2); ++i) { //back
		for (int j = 0; j < ((mp)+5); ++j) {
			fprintf(writeGrid, "2 %lf %lf %lf\n", ((i*dx) + 2 * dx) + offsetX, (2 * dy + j*dy) + offsetY, (op + 1)*dz + offsetZ);
		}
	}

	//PARTICULAS TIPO 3 #interna
	for (int k = 0; k < (op + 6); k++) {
		for (int i = 0; i < (2 * mp) + 6; ++i) { //down
			fprintf(writeGrid, "3 %lf %lf %lf\n", (i*dx) + offsetX, offsetY, dz*k - 3 * dz + offsetZ);
		}
		for (int i = 0; i < (mp)+6; ++i) { //left
			fprintf(writeGrid, "3 %lf %lf %lf\n", offsetX, ((i*dy) + dy) + offsetY, dz*k - 3 * dz + offsetZ);
		}
		for (int i = 0; i < (mp)+6; ++i) { //right
			fprintf(writeGrid, "3 %lf %lf %lf\n", (2 * xl + (5 * dx)) + offsetX, (i*dy + dy) + offsetY, dz*k - 3 * dz + offsetZ);
		}
	}

	for (int i = 0; i < ((2 * mp) + 4); ++i) { //front
		for (int j = 0; j < ((mp)+6); ++j) {
			fprintf(writeGrid, "3 %lf %lf %lf\n", i*dx + dx + offsetX, j*dy + dy + offsetY, -3 * dz + offsetZ);
		}
	}
	for (int i = 0; i < ((2 * mp) + 4); ++i) { //back
		for (int j = 0; j < ((mp)+6); ++j) {
			fprintf(writeGrid, "3 %lf %lf %lf\n", i*dx + dx + offsetX, j*dy + dy + offsetY, (op + 2)*dz + offsetZ);
		}
	}


	fclose(writeGrid);
}


int main()
{
	int nump = 0;

	data_in *input_data = new data_in();
	//Particle2D *particles;
	Particle3D *particles;
	ReadWrite FileControl;

	//criaParticulasQuebraBarragem3D();
	criaParticulasQuebraBarragem3D_offset_semzero();


	// Pegando e imprimindo número de particulas "nump"
	ifstream in;
	in.open("mps.grd");
	while (!in.is_open()){
		std::cout << "Tentando ler mps.grd novamente" << endl;
		in.open("mps.grd");
	}
	in >> nump; std::cout << nump << endl;

	// Array de partículas
	//particles = new Particle2D[nump];
	particles = new Particle3D[nump];

	// Lendo informações sobre as partículas e a simulação
	FileControl.ReadData("mps.data", input_data);
	FileControl.ReadParticles3D(particles, "mps.grd", nump);

	// Pegando o valor do passo de tempo
	input_data->dt = input_data->MaxDt;

	// Guardando no vtu0 a configuação inicial do sistema/simulação
	FileControl.WriteOut(particles, 0, nump);

	TimeIntegration(particles, FileControl, input_data, nump);

	delete[] particles;
	return 0;
}