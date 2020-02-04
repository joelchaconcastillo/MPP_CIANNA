//Estructura que almacena informacion de los platos 
struct infoPlatos {
	string nombre;				  //Nombre del plato
	float precio;           //Precio del plato
	float cantidad;					//Cantidad en gramos
	vector<float> infoN;    //Informacion nutricional del plato
	vector<string> alg;     //Alergenos del plato
	vector<string> inc;     //Incompatibilidades alimenticias del plato
	vector<int> gruposAl;
};

