#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <armadillo>
#include <string>
#include <math.h>

using namespace std;
using namespace arma;


struct Shells{
    double exp;
    double coord_1D;
    vector<double> coord;
    int total_L; //total angular momentum numbers
    mat A; //angular quantum number matrix
    double stdev;
};


//set angular momentum matrix will all possible orientations given total quantum number
void cartesian_angular_momentum(mat &A, const int &total_L){
    int total_sum;

    for(int i = 0; i <= total_L; i++){
        for(int j = 0; j <= total_L; j++){
            for(int k = 0; k <= total_L; k++){
                total_sum = i + j + k;
                if (total_sum == total_L){
                    A.resize(A.n_rows + 1, 3);
                    A(A.n_rows - 1, 0) = k;
                    A(A.n_rows - 1, 1) = j;
                    A(A.n_rows - 1, 2) = i;
                }
            }
        }
    }
}


//counts number of elements from input line
int element_counter(const string line){
    istringstream isstring(line);
    double temp;
    int count = 0;

    while(isstring >> temp){
        count++;
    }
    return count;
}


//reads file and sets value of shells in a vector of shells
void read_shell(vector<Shells>& Shell, const string& file_name){
    ifstream in;

    cout << '\n';
    in.open(file_name, ios::in); //opening file + catch errors
    if(!in.is_open()){
        throw runtime_error(file_name + " didn't open!");
    } else {
        cout << "file opened!" << endl;
    }

    string line;
    while (getline(in, line)) { //iteration through every line in file
        int num_of_elements = element_counter(line);
        istringstream isstring(line);
        Shells input_shell;
        string number_str;

        if (num_of_elements == 3){ //1D coordinate input read
            isstring >> input_shell.coord_1D >> input_shell.exp >> input_shell.total_L;

            isstring >> input_shell.coord_1D;
            isstring >>input_shell.exp;
            isstring >> input_shell.total_L;

        } else if (num_of_elements > 3){ //larger coordinate dim input read
            for(int i = 0; i < num_of_elements - 2; i++){
                isstring >> number_str;
                input_shell.coord.push_back(stod(number_str));
            }
            isstring >> input_shell.exp >> input_shell.total_L;
        } else{ //less than 1D coordinate input (file format error)
            throw invalid_argument("Problem with data format");
        }
        cartesian_angular_momentum(input_shell.A, input_shell.total_L); //set angular momentum orientations
        input_shell.stdev = pow(sqrt(2.0 * input_shell.exp), -1.0); //set standard dev of gaussian
        Shell.push_back(input_shell);
    }
    in.close();
}


//evaluates and returns Gaussian function
double g(double x, Shells shell){
    return pow((x - shell.coord_1D), shell.total_L) * exp(-1 * shell.exp * pow((x - shell.coord_1D), 2));
}


//returns 1D overlap of two gaussian integrals
double overlap_1D(double x, Shells shell_a, Shells shell_b){
    return g(x, shell_a) * g(x, shell_b);
}


//returns evaluation of trapazoidal integral with adjusted bounds
double ext_trap_num_int(Shells shell_a, Shells shell_b, int n, double epsilon, double correction_coeff){
    double upper, lower, upper_mag, lower_mag, I, I_previous, f_x_0, f_x_n, f_x_i, h;

    // initial I calculation
    if(shell_a.coord_1D - shell_a.stdev <= shell_b.coord_1D - shell_b.stdev){
        lower_mag = shell_a.stdev;
        lower = shell_a.coord_1D - lower_mag;
    } else {
        lower_mag = shell_b.stdev;
        lower = shell_b.coord_1D - lower_mag;
    }

    if(shell_a.coord_1D + shell_a.stdev <= shell_b.coord_1D + shell_b.stdev){
        upper_mag = shell_b.stdev;
        upper = shell_b.coord_1D + upper_mag;
    } else {
        upper_mag = shell_a.stdev;
        upper = shell_a.coord_1D + upper_mag;
    }

    I_previous = 0.0;
    h = (upper - lower) / n;
    I = 0.0;
    f_x_0 = overlap_1D(lower, shell_a, shell_b);
    f_x_n = overlap_1D(upper, shell_a, shell_b);

    for (int i = 1; i < n; i++) {
        f_x_i = overlap_1D(lower + (h * i), shell_a, shell_b);
        I += h * (f_x_0 + f_x_n + (2 * f_x_i)) / 2;
    }


    //loop to find bound
    int j = 0;
    while(!(abs(I - I_previous) < epsilon)){
        //cout << "Iteration: " << j << '\n';
        I_previous = I;

        lower = lower - (correction_coeff * lower_mag);
        upper = upper + (correction_coeff * upper_mag);

        h = (upper - lower) / n;
        I = 0.0;
        f_x_0 = overlap_1D(lower, shell_a, shell_b);
        f_x_n = overlap_1D(upper, shell_a, shell_b);

        for (int i = 1; i < n; i++) {
            f_x_i = overlap_1D(lower + (h * i), shell_a, shell_b);
            I += h * (f_x_0 + f_x_n + (2 * f_x_i)) / 2;
        }
        j++;
    }

    return I;
}


int factorial(int l){
    int total = 1;

    for(int i = l; i > 0; i--){
        total *= i;
    }
    return total;
}

int double_factorial(int l){
    int total = 1;

    if((l % 2) == 0){
        for(int i = l; i > 0; i--){
            if((i % 2) == 0){
                total *= i;
            }
        }
    } else{
        for(int i = l; i > 0; i--){
            if((i % 2) == 1){
                total *= i;
            }
        }
    }
    return total;
}

int binomial(int m, int n){
    return factorial(m) / (factorial(n) * factorial(m - n));
}

//returns overlap integration with analytical integration in 1D
double overlap_1D_ana_int(double exp_a, double exp_b, int l_a, int l_b, double coord_a, double coord_b){
    double sum = 0.0;
    double coord_p = ((exp_a * coord_a) + (exp_b * coord_b)) / (exp_a + exp_b);
    double prefactor = exp((-1 * exp_a * exp_b * pow((coord_a - coord_b),2)) / (exp_a + exp_b));

    for(int i = 0; i <= l_a; i++){
        for(int j = 0; j <= l_b; j++){
            if(((i + j) % 2) == 0){
                sum += binomial(l_a, i) * binomial(l_b, j) * ((double_factorial(i + j - 1) *
                                                               pow((coord_p - coord_a),(l_a - i)) * pow((coord_p - coord_b),(l_b - j))) /
                                                              pow((2 * (exp_a + exp_b)),((i + j) * 0.5)));
            }
        }
    }
    return prefactor * sqrt((M_PI) / (exp_a + exp_b)) * sum;
}


//returns overlap integration with analytical in 3D
double overlap_3D_ana_int(double exp_a, double exp_b, rowvec l_a, rowvec l_b, vector<double> coord_a, vector<double> coord_b){
    double total_overlap = 1;

    for(int i = 0; i < 3; i++){
        total_overlap *= overlap_1D_ana_int(exp_a, exp_b, l_a[i], l_b[i], coord_a[i], coord_b[i]);
    }
    return total_overlap;
}


//returns matrix of overlap integral values
mat analytical_overlap_int(Shells shell_a, Shells shell_b){
    mat overlap_integral(shell_a.A.n_rows, shell_b.A.n_rows);

    for(int i = 0; i < shell_a.A.n_rows; i++){
        for(int j = 0; j < shell_b.A.n_rows; j++){
            overlap_integral(i, j) = overlap_3D_ana_int(shell_a.exp, shell_b.exp, shell_a.A.row(i), shell_b.A.row(j), shell_a.coord, shell_b.coord);
        }
    }

    return overlap_integral;
}


int main() {

    //input file directory
    vector<string> numerical_input_fname{"/Users/joey/Desktop/Chem_179/homework_2/sample_input/numerical/1.txt",
                                   "/Users/joey/Desktop/Chem_179/homework_2/sample_input/numerical/2.txt",
                                   "/Users/joey/Desktop/Chem_179/homework_2/sample_input/numerical/3.txt"};

    vector<string> analytical_input_fname{"/Users/joey/Desktop/Chem_179/homework_2/sample_input/analytical/1.txt",
                                   "/Users/joey/Desktop/Chem_179/homework_2/sample_input/analytical/2.txt",
                                   "/Users/joey/Desktop/Chem_179/homework_2/sample_input/analytical/3.txt"};

    //output for question 1
    for(int i = 0; i < 3; i++){
        vector<Shells> test;
        read_shell(test, numerical_input_fname[i]);
        cout << "Numerical Integration " << i + 1 << ".txt" << '\n';

        // numerical integral stepsize
        int n = 1000;

        cout << "Shells A: " << '\n';
        cout << "coord_1D: " << test[0].coord_1D << '\n';
        cout << "exp: " << test[0].exp << '\n';
        cout << "quantum numbers " << '\n';
        cout << test[0].A << '\n';

        cout << "Shells B: " << '\n';
        cout << "coord_1D: " << test[1].coord_1D << '\n';
        cout << "exp: " << test[1].exp << '\n';
        cout << "quantum numbers " << '\n';
        cout << test[1].A << '\n';

        cout << "ID numerical Integral: " << ext_trap_num_int(test[0], test[1], n,
                                                              pow(10, -7), 1.5);
        cout << '\n' << '\n';
    }


    //output for question 2
    for(int i = 0; i < 3; i++){
        vector<Shells> test;
        read_shell(test, analytical_input_fname[i]);
        cout << "Analytical Integration " << i + 1 << ".txt" << '\n';

        cout << "Shells A: " << '\n';
        cout << "coord: ";
        for(int j = 0; j < 3; j++){
            cout << test[0].coord[j] << ' ';
        }
        cout << '\n';
        cout << "exp: " << test[0].exp << '\n';
        cout << "total quantum numer: " << test[0].total_L << '\n';
        cout << "quantum numbers " << '\n';
        cout << test[0].A << '\n';

        cout << "Shells B: " << '\n';
        cout << "coord: ";
        for(int j = 0; j < 3; j++){
            cout << test[1].coord[j] << ' ';
        }
        cout << '\n';
        cout << "exp: " << test[1].exp << '\n';
        cout << "total quantum numer: " << test[1].total_L << '\n';
        cout << "quantum numbers " << '\n';
        cout << test[1].A << '\n';
        cout << '\n';

        cout << "3D Overlap Analytical Integral: " << '\n';
        cout << '\n' << analytical_overlap_int(test[0], test[1]) << '\n';
    }

    cout << '\n' << "***Order of Analytical Overlap Integral corresponds to order of quantum numbers listed in each case.***" << '\n';

    return 0;
}

