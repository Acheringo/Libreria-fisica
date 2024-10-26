#ifndef fisicaUZ_H_INCLUDED
#define fisicaUZ_H_INCLUDED

#include<iostream>

using namespace std;

class complex{ //Clase de números complejos
    public:
        double re; //Parte real
        double im; //Parte imaginaria
        void print(){ //Método para mostrar en pantalla el número complejo, formateado como re+im i
            if(re==0 && im==0){
                cout<<0;
            }
            if(re<0){
                cout<<re;
            }
            if(re>0){
                cout<<re;
            } // Solo mostramos la parte real si es distinta de 0
            if(re==0){
                if(im>0){
                    cout<<im<<"i";
                }
                if(im<0){
                    cout<<im<<"i";
                }
                //Solo mostramos la parte imaginaria si es distinta de 0
            }
            if(re!=0){
                if(im>0){
                    if(im==1){
                        cout<<"+i";
                    }
                    else{
                        cout<<"+"<<im<<"i";
                    }
                } //Si la parte imaginaria es mayor que 0, la mostramos con un +
                if(im<0){
                    if(im==-1){
                        cout<<"-i";
                    }
                    else{
                    cout<<im<<"i";
                    }
                } //Si la parte imaginaria es menor que 0, la mostramos sola, ya que el signo ya está puesto
            }
        }
        complex(){ // Constructor por defecto
            re = 1.0;
            im = 0.0;
        }
        complex(double a, double b){ // Constructor en forma binómica
            re = a;
            im = b;
        }
        complex(){ \\ Constructor por defecto
            re = 1.0;
            im = 0.0;
        }
        complex(double a, double b){ \\ Constructor en forma binómica
            re = a;
            im = b;
        }
        double mod(){ //Método para calcular el cuadrado del módulo de un número complejo
            return re*re+im*im;
        }
        double tg(){ //Método que devuelve el arcotagente del ángulo polar del número complejo
            return im/re;
        }
        complex operator*(const complex& other) {
            complex result;
            result.re = (re * other.re) - (im * other.im);
            result.im = (re * other.im) + (im * other.re);
            return result;
        }
        complex operator+(const complex& other) {
            complex result;
            result.re = re+other.re;
            result.im = im+other.im;
            return result;
        }
        complex operator/(const complex& other){
            complex res;
            double denom=other.re*other.re+other.im*other.im;
            res.re=(re*other.re+im*other.im)/denom;
            res.im=(im*other.re-re*other.im)/denom;
            return res;
        }
        complex& operator=(const complex& other){
            re=other.re;
            im=other.im;
            return *this;
        }

};

class matrix{ //Clase de matrices reales
    public:
        int nrow; //Número de filas
        int ncol; //Número de columnas
        double** entrie; //Entradas, puntero doble (fila y columna)

        matrix(int rows, int cols){ //Constructor dinámico
            nrow = rows; //El número de filas que se introduce es el de la matriz
            ncol = cols; //El número de columnas que se introduce es el de la matriz
            entrie = new double*[nrow]; //Creamos un array de punteros, con un número de elementos igual al número de filas
            for (int i = 0; i < nrow; ++i){
                entrie[i] = new double[ncol]; //Para cada fila, creamos un array dentro del que ya habíamos creado antes, tantos como sea el número de columnas
            }
        }

        ~matrix(){ //Destructor de la matriz
            for (int i = 0; i < nrow; ++i){
                delete[] entrie[i]; //Borramos las columnas
            }
            delete[] entrie; //Borramos las filas
        }
            
        matrix(const matrix& other){ // Constructor de copia (copia profunda)
            nrow = other.nrow; //Asignamos el mismo número de filas
            ncol = other.ncol; //Asignamos el mismo número de columnas
            entrie = new double*[nrow]; //Creamos filas
            for (int i = 0; i < nrow; ++i){
                entrie[i] = new double[ncol]; //Creamos columnas
                for (int j = 0; j < ncol; ++j){
                    entrie[i][j] = other.entrie[i][j]; //Copiamos elemento a elemento
                }
            }
        }


        matrix& operator=(const matrix& other){ //Sobrecarga del operador de asignación (también copia profunda)
            if (this != &other) {
                //Liberar la memoria actual para crear una nueva
                for (int i = 0; i < nrow; ++i){
                    delete[] entrie[i];
                }
                delete[] entrie;

                //Proceso de copia anterior
                nrow = other.nrow;
                ncol = other.ncol;
                entrie = new double*[nrow];
                for (int i = 0; i < nrow; ++i){
                    entrie[i] = new double[ncol];
                    for (int j = 0; j < ncol; ++j){
                        entrie[i][j] = other.entrie[i][j];
                    }
                }
            }
            return *this;
        }

        matrix operator*(const matrix& other){ //Sobrecarga del operador de producto
            matrix res(nrow, other.ncol);
            if(ncol!=other.nrow){
                cout<<"Matrices dimensions are incompatible";
            }
            else{
                for (int i = 0; i < res.nrow; i++){
                    for (int j = 0; j < res.ncol; j++){
                        res.entrie[i][j] = 0;
                        for (int r = 0; r < ncol; r++){
                            res.entrie[i][j] += entrie[i][r] * other.entrie[r][j]; //Hacemos la multiplicación elemento por elemento (siguiendo las normas)
                        }
                    }
                }
            }
            return res;
        }

        matrix operator+(const matrix& other){ //Sobrecarga del operador de suma
            matrix res(nrow, ncol);
            if(ncol!=other.ncol || nrow!=other.nrow){
                cout<<"Matrices dimensions are incompatible";
            }
            else{
                for (int i = 0; i < nrow; i++){
                    for (int j = 0; j < res.ncol; j++){
                        res.entrie[i][j]=entrie[i][j]+other.entrie[i][j];
                    }
                }
            }
            return res;
        }
        void print(){ //Método para mostrar los valores en pantalla
            for(int i=0; i<nrow; i++){
                for(int j=0; j<ncol;j++){
                    cout<<entrie[i][j]<<" ";
                }
                cout<<endl;
            }
        }
};

class complexmatrix{ //Clase de matrices reales (iguales que matrix solo que con entradas complex en lugar de double)
    public:
        int nrow;
        int ncol;
        complex** entrie;

        complexmatrix(int rows, int cols){
            nrow = rows;
            ncol = cols;
            entrie = new complex*[nrow];
            for (int i = 0; i < nrow; ++i){
                entrie[i] = new complex[ncol];
            }
        }

        ~complexmatrix(){
            for(int i=0; i<nrow; i++){
                delete[] entrie[i];
            }
            delete[] entrie;
        }

        complexmatrix(const complexmatrix& other){
            nrow = other.nrow;
            ncol = other.ncol;
            entrie = new complex*[nrow];
            for (int i = 0; i < nrow; ++i){
                entrie[i] = new complex[ncol];
                for (int j = 0; j < ncol; ++j){
                    entrie[i][j] = other.entrie[i][j];
                }
            }
        }

        complexmatrix& operator=(const complexmatrix& other){
            if (this != &other) {
                for (int i = 0; i < nrow; ++i){
                    delete[] entrie[i];
                }
                delete[] entrie;

                nrow = other.nrow;
                ncol = other.ncol;
                entrie = new complex*[nrow];
                for (int i = 0; i < nrow; ++i){
                    entrie[i] = new complex[ncol];
                    for (int j = 0; j < ncol; ++j){
                        entrie[i][j] = other.entrie[i][j];
                    }
                }
            }
            return *this;
        }
        complexmatrix operator*(const complexmatrix& other){
            complexmatrix res(nrow, other.ncol);
            if(ncol!=other.nrow){
                cout<<"Matrices dimensions are incompatible";
            }
            else{
                for (int i = 0; i < res.nrow; i++){
                    for (int j = 0; j < res.ncol; j++){
                        res.entrie[i][j].re = 0;
                        res.entrie[i][j].im=0;
                        for (int r = 0; r < ncol; r++){
                            res.entrie[i][j] = res.entrie[i][j] + entrie[i][r]*other.entrie[r][j]; //Hacemos la multiplicación elemento por elemento (siguiendo las normas)
                        }
                    }
                }
            }
            return res;
        }

        complexmatrix operator+(const complexmatrix& other){
            complexmatrix res(nrow, ncol);
            if(ncol!=other.ncol || nrow!=other.nrow){
                cout<<"Matrices dimensions are incompatible";
            }
            else{
                for (int i = 0; i < nrow; i++){
                    for (int j = 0; j < res.ncol; j++){
                        res.entrie[i][j]=entrie[i][j]+other.entrie[i][j];
                    }
                }
            }
            return res;
        }
        void print(){
            for(int i=0; i<nrow; i++){
                for(int j=0; j<ncol;j++){
                    entrie[i][j].print(); cout<<" ";
                }
                cout<<endl;
            }
        }
};

class tensor { //Clase de tensores
private:
    double* data; //Datos
    int* dimensions; //Dimensiones del tensor
    int order; //Orden del tensor
    int totalSize; //Tamaño total del tensor

    int calculateIndex(int* indices) const{ // Método para calcular el índice lineal (la explicación es larga, está al finalizar el método)
        int index = 0;
        int multiplier = 1;
        for (int i = order - 1; i >= 0; --i){
            index += indices[i] * multiplier;
            multiplier *= dimensions[i];
        }
        return index;
    }
    /*
    Vamos a estar guardando los elementos del tensor en un array lineal, ya que c++ no tiene soporte para arrays de más de 2 dimensiones
    Para ello, vamos a estar recorriendo los índices desde el final hasta el comienzo. A cada elemento del tensor le corresponderá
    un elemento en el array que será la suma de los índices por un multiplicador, que este será el producto acumulado de las dimensiones,
    de esta forma, hacemos una asignación 1 a 1 de un elemento del tensor con un elemento del array.
    */


    void printRecursive(int* indices, int depth) const{//Método para mostrar el tensor (explicación larga)
        if (depth == order){
            int index = calculateIndex(indices);
            cout << data[index] << " ";
            return;
        }

        cout << "[";
        for (int i = 0; i < dimensions[depth]; ++i){
            indices[depth] = i;
            printRecursive(indices, depth + 1);
            if (i < dimensions[depth] - 1){
                cout << ", ";
            }
        }
        cout << "]"<<endl;
        if (depth == 0){
           cout << endl;
        }
    }
    /*
    Este método es recursivo, es decir, se llama a sí mismo. Al introducir unos índices y una profundidad
    (algo así como la capa de matrices), va mostrando en pantalla los elementos del tensor, formateados
    en forma de hipermatriz para facilitar la lectura. Hace uso de calculateIndex para saber qué elemento
    tiene que mostrar cada vez, y se llama así mismo hasta que llegue a la última capa (que es igual que order)
    */
    public:
        tensor(int* dims, int ord) : order(ord){ //Constructor dinámico
            dimensions = new int[order]; //Las dimensiones del tensor van a ser un array con "order" elementos
            totalSize = 1; //Inicializamos el tamaño total como 1
            for (int i = 0; i < order; ++i) {
                dimensions[i] = dims[i]; //A cada elemento de dimensiones, le asignamos el valor del array de dimensiones que entra
                totalSize *= dimensions[i]; //El tamaño total quedará multiplicado por cada una de las dimensiones
            }
            data = new double[totalSize]; //El array con los datos tendrá un número de elementos igual a tamaño total
        }

        // Destructor
        ~tensor(){
            delete[] data; //Borra los datos
            delete[] dimensions; //Borra las dimensiones
        }

        tensor(const tensor& other){
            order = other.order;
            dimensions = new int[order];
            totalSize = other.totalSize;
            for (int i = 0; i < order; ++i){
                dimensions[i]=other.dimensions[i];
            }
            data=new double[totalSize];
            for(int i=0; i<totalSize; i++){
                data[i]=other.data[i];
            }
        }

        tensor& operator=(const tensor& other){
            if (this != &other) {
                delete[] data;
            
                order = other.order;
                dimensions=new int[order];
                totalSize = other.totalSize;
                data = new double[totalSize];
                for(int i=0; i<order;i++){
                    dimensions[i]=other.dimensions[i];
                }
                for (int i = 0; i < totalSize; ++i){
                    data[i]=other.data[i];
                }
            }
            return*this;
        }
        void setValue(int* indices, double value){ //Méotodo para asignar valores
            if (indices != nullptr){ //Si el array de índices no es nulo 
                int index = calculateIndex(indices); //Calculamos la posición del array de datos
                data[index] = value; //Asignamos el valor a esa posición
            }
        }

        double getValue(int* indices) const{ //Método para obtener un valor
        try{
            if (indices == nullptr){ //Si el array de índices es nulo
                throw invalid_argument("Element does not exist"); //Mostramos mensaje de error
            }
            int index = calculateIndex(indices); //Calculamos la posición en el array de datos
            return data[index]; //Devolvemos el elemento
        } catch(const invalid_argument&e){ //Error
            cout<<e.what()<<endl;
            }
            return 0;
        }
        
        void print() const{ //Método para imprimir el tensor
            int* indices = new int[order]; //Creamos un array para los índices del tamaño del orden del tensor
            printRecursive(indices, 0); //Iniciamos el print recursivo
            delete[] indices; //Borramos los índices
        }
        void printVal(int* indices) const{ //Método para imprimir un elemento
            cout<<getValue(indices); //Muestra el elemento correspondiente a los índices introducidos
        }
        tensor operator*(const tensor& other) const {
            // El orden del nuevo tensor es la suma de los órdenes de ambos tensores
            int newOrder = this->order + other.order;
            
            // Crear un array con las nuevas dimensiones combinadas
            int* newDimensions = new int[newOrder];
            for (int i = 0; i < this->order; ++i) {
                newDimensions[i] = this->dimensions[i];
            }
            for (int i = 0; i < other.order; ++i) {
                newDimensions[this->order + i] = other.dimensions[i];
            }
            
            // Crear el nuevo tensor
            tensor result(newDimensions, newOrder);
            delete[] newDimensions;
            
            // Rellenar los datos del nuevo tensor aplicando el producto tensorial
            // Recorremos ambos tensores y asignamos los valores al nuevo tensor
            int* indicesA = new int[this->order];
            int* indicesB = new int[other.order];
            int* indicesResult = new int[newOrder];
            
            for (int i = 0; i < this->totalSize; ++i) {
                // Obtener los índices correspondientes del tensor A
                int remaining = i;
                for (int j = this->order - 1; j >= 0; --j) {
                    indicesA[j] = remaining % this->dimensions[j];
                    remaining /= this->dimensions[j];
                }
                
                for (int j = 0; j < other.totalSize; ++j) {
                    // Obtener los índices correspondientes del tensor B
                    remaining = j;
                    for (int k = other.order - 1; k >= 0; --k) {
                        indicesB[k] = remaining % other.dimensions[k];
                        remaining /= other.dimensions[k];
                    }
                    
                    // Combinar los índices de A y B en los índices del nuevo tensor
                    for (int a = 0; a < this->order; ++a) {
                        indicesResult[a] = indicesA[a];
                    }
                    for (int b = 0; b < other.order; ++b) {
                        indicesResult[this->order + b] = indicesB[b];
                    }
                    
                    // Asignar el producto de los valores de A y B al nuevo tensor
                    double value = this->getValue(indicesA) * other.getValue(indicesB);
                    result.setValue(indicesResult, value);
                }
            }
            
            delete[] indicesA;
            delete[] indicesB;
            delete[] indicesResult;
            
            return result;
        }

};

class poly {
    public:
        int dim;
        double* entrie;
        poly(int dimension){
            dim=dimension;
            entrie=new double[dim];
        }
        ~poly(){
            delete[] entrie;
        }
        poly(const poly& other){
            dim = other.dim;
            entrie = new double[dim];
            for (int i = 0; i < dim; ++i){
                entrie[i]=other.entrie[i];
            }
        }

        poly& operator=(const poly& other){
            if (this != &other) {
                delete[] entrie;
            
                dim = other.dim;
                entrie = new double[dim];
                for(int i=0; i<dim;i++){
                    entrie[i]=other.entrie[i];
                }
            }
            return*this;
        }
        poly operator+(const poly&other){
            int dmax;
            int dmin;
            if(dim>other.dim){
                dmax=dim;
                dmin=other.dim;
                poly res(dmax);
                for(int i=0; i<dmin;i++){
                    res.entrie[i]=entrie[i]+other.entrie[i];
                }
                for(int i=0; i<dmax-dmin;i++){
                    res.entrie[i+dmin]=entrie[i+dmin]+0;
                }
                return res;
            }
            else{
                dmax=other.dim;
                dmin=dim;
                poly res(dmax);
                for(int i=0; i<dmin;i++){
                    res.entrie[i]=entrie[i]+other.entrie[i];
                }
                for(int i=0; i<dmax-dmin;i++){
                    res.entrie[i+dmin]=other.entrie[i+dmin]+0;
                }
                return res;
            }
        }

        void print(){
            for(int i=0; i<dim; i++){
                if(i==0){
                    cout<<entrie[i];
                }
                if(i==1){
                    if(entrie[i]>0){
                        cout<<"+"<<entrie[i]<<"x";
                    }
                    if(entrie[i]<0){
                        cout<<entrie[i]<<"x";
                    }
                }
                if(i>1){
                    if(entrie[i]>0){
                        cout<<"+"<<entrie[i]<<"x^"<<i;
                    }
                    if(entrie[i]<0){
                        cout<<entrie[i]<<"x^"<<i;
                    }
                }
            }
        }
        double evaluate(int x){
            double res=0;
            for(int i=0; i<dim; i++){
                double pow=1;
                for(int j=1; j<i; j++){
                    pow*=x;
                }
                if(x==0){
                    res=entrie[0];
                }
                else{
                    res+=entrie[i]*pow;
                }
            }
            return res;
        }
};

complex conjugate(complex a){ //Complejo conjugado
    /*
    conjugado de a+bi=a-bi
    */
    complex res;
    res.re=a.re;
    res.im=-a.im;
    return res;
}

int kronecker(int i, int j){ //Delta de kronecker
    if(i==j){ //Si i es igual que j, la delta vale 1
        return 1;
    }
    else{ //Si son distintos, vale 0
        return 0;
    }
}

int levi_civita(int i, int j, int k){ //Tensor de Levi-Civita
    if(i==j || i==k || j==k){ //Si hay dos índices repetidos, vale 0
        return 0;
    }
    else{ //Explicación larga
        int perm=1;
        int ijk[3]={i,j,k};
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 3-i-1; j++) {
                if (ijk[j] > ijk[j+1]) {
                    int temp = ijk[j];
                    ijk[j] = ijk[j+1];
                    ijk[j+1] = temp;
                    perm*=(-1);
                }
            }
        }
        return perm;
    }
    /*
    El tensor vale 1 si los índices están en orden cíclico, y -1 si es anticíclico. Para calcularlo vamos a ver el número
    de permutaciones que tenemos que hacerle a los índices para que estén ordenados, ya que la definición que hemos dado
    es equivalente a la signatura de las permutaciones que se le ha hecho a los índices desde la posición ordenada (de
    menor a mayor). Por tanto, por cada movimiento que hagamos hasta que estén ordenados, multiplicamos por -1 el valor
    de la signatura (llamado perm), y cuando esté ordenado, lo devolvemos
    */
}

int factorial(int n){ //Número factorial (no recursivo)
    int resultado=1;
    for(int i=1; i<=n; i++){
        resultado*=i; //Multiplica por todos los números desde 1 hasta el introducido (definición de factorial)
    }
    return resultado;
}

int comb(int n, int k){ //Número combinatorio n sobre k
    int resultado;
    resultado = (factorial(n)/(factorial(k)*factorial(n-k))); //Definición de número combinatorio
    return resultado;
}

double prod_esc(double arr1[], double arr2[], int n){ //Producto escalar de dos vectores
    double res = 0;
    for(int i = 0; i < n; i++){
        res += arr1[i] * arr2[i]; //Definición del producto escalar de dos vectores en el espacio euclídeo
    }
    return res;
}

void iseq(int n, int m){ //Función que devuelve error si dos números son distintos, usado en matrices
    if (n != m) {
        throw invalid_argument("Matrix dimension incompatible");
    }
}

matrix gaussJordan(matrix&mat){ //Función que trianguliza la matriz por el método de Gauss-Jordan 
    matrix result = mat; //Copiar la matriz original

    for (int i = 0; i < result.nrow; i++){
        if (i >= result.ncol) break;  //Si la matriz no es cuadrada, evitamos el desbordamiento

        //Buscamos el máximo de la columna i
        double maxEl = result.entrie[i][i];
        int maxRow = i;

        //Busca números mayores en las filas de debajo
        for (int k = i + 1; k < result.nrow; k++){
            if (result.entrie[k][i] > maxEl){
                maxEl = result.entrie[k][i];
                maxRow = k;
            }
        }

        //Intercambia la fila si es necesario
        if (i != maxRow){
            for (int k = 0; k < result.ncol; k++){
                swap(result.entrie[maxRow][k], result.entrie[i][k]);
            }
        }

        if (result.entrie[i][i] == 0) {
            return result;  //Maneja el caso de que tengamos ceros en la diagonal
        }

        //Elimina elementos debajo del pivote
        for (int k = i + 1; k < result.nrow; k++) {
            double factor = result.entrie[k][i] / result.entrie[i][i];
            for (int j = i; j < result.ncol; j++) {
                result.entrie[k][j] -= factor * result.entrie[i][j];
            }
        }
    }
    return result; //Devolvemos la matriz triangulada
}

complexmatrix cgaussJordan(complexmatrix& mat){ //Lo mismo con matrices complejas
    complexmatrix result = mat;

    for (int i = 0; i < result.nrow; ++i) {
        double maxMod = result.entrie[i][i].mod();
        int maxRow = i;

        for (int k = i + 1; k < result.nrow; ++k) {
            if (result.entrie[k][i].mod() > maxMod) {
                maxMod = result.entrie[k][i].mod();
                maxRow = k;
            }
        }

        if (i != maxRow) {
            for (int k = 0; k < result.ncol; ++k) {
                swap(result.entrie[maxRow][k], result.entrie[i][k]);
            }
        }

        if (result.entrie[i][i].mod() == 0) {
            return result;
        }

        for (int k = i + 1; k < result.nrow; ++k) {
            complex factor = result.entrie[k][i]/result.entrie[i][i];

            for (int j = i; j < result.ncol; ++j) {
                complex prodFactor = factor* result.entrie[i][j];
                result.entrie[k][j].re -= prodFactor.re;
                result.entrie[k][j].im -= prodFactor.im;
            }
        }
    }

    return result;
}

double det(matrix mat){ //Función que calcula el determinante de una matriz
    matrix gauss=gaussJordan(mat); //Obtenemos la matriz triangulada de la introducida
    double res=1;
    try{
        iseq(mat.ncol, mat.nrow); //Miramos si la matriz es cuadrada
        for(int i=0; i<mat.ncol; i++){
            res*=gauss.entrie[i][i]; //Multiplicamos el resultado por los elementos de la diagonal
        }
    } catch (const invalid_argument& e){
        cout << e.what() << endl; //Error
    }
    return res;
}

complex cdet(complexmatrix mat){ //Lo mismo para matrices complejas
    complexmatrix gauss=cgaussJordan(mat);
    complex res;
    res.im=1;
    res.re=1;
    try{
        iseq(mat.ncol, mat.nrow);
        for(int i=0; i<mat.ncol; i++){
            res=res*gauss.entrie[i][i];
        }
    } catch (const invalid_argument& e) {
        cout << e.what() << endl;
    }
    return res;
}

poly derivative(poly p){
    poly res(p.dim-1);
    for(int i=1; i<p.dim; i++){
        res.entrie[i-1]=p.entrie[i]*i;
    }
    return res;
}
#endif //fisicaUZ_H_INCLUDED