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
        Shell.push_back(input_shell);
    }
    in.close();
}


//returns function evaluation of overlap guassian
double gauss_overlap_1D(const double x_a, const double x_b, const double alpha, const double beta, const int l_a,
                        const int l_b, const double x){
    return pow((x - x_a), l_a) * pow((x - x_b), l_b) *
    exp((-1 * alpha * pow((x - x_a), 2)) + (-1 * beta * pow((x - x_b), 2)));
}



double ext_trapazoidal_int(const double x_a, const double x_b, const double alpha, const double beta, const int l_a,
                           const int l_b, const double x, double (*func)(double, double, double, double, int, int, double)){
    double std_A = pow(sqrt(2.0 * alpha), -1.0);
    double std_B = pow(sqrt(2.0 * beta), -1.0);
    double lower = 4.0 * std_A
}




int main() {

    vector<string> numerical_fname{"/Users/joey/Desktop/Chem_179/homework_2/sample_input/numerical/1.txt",
                                   "/Users/joey/Desktop/Chem_179/homework_2/sample_input/numerical/2.txt",
                                   "/Users/joey/Desktop/Chem_179/homework_2/sample_input/numerical/3.txt"};

    vector<string> analytical_fname{"/Users/joey/Desktop/Chem_179/homework_2/sample_input/analytical/1.txt",
                                   "/Users/joey/Desktop/Chem_179/homework_2/sample_input/analytical/2.txt",
                                   "/Users/joey/Desktop/Chem_179/homework_2/sample_input/analytical/3.txt"};

    vector<Shells> test_1;
    read_shell(test_1, numerical_fname[0]);


    cout << '\n' << '\n' << "1.txt Shells A: " << '\n';
    cout << "coord_1D: " << test_1[0].coord_1D << '\n';
    cout << "coord: " << '\n';
    for(int i = 0; i < test_1[0].coord.size(); i++){
        cout << test_1[0].coord[i] << ' ';
    }
    cout << '\n' << "exp: " << test_1[0].exp << '\n';
    std::cout << "Data type of exp: " << typeid(test_1[0].exp).name() << std::endl;
    cout << "quantum numbers " << '\n';
    cout << test_1[0].A << '\n';



    cout << '\n' << "1.txt Shells B: " << '\n';
    cout << "coord_1D: " << test_1[1].coord_1D << '\n';
    cout << "coord: " << '\n';
    for(int i = 0; i < test_1[0].coord.size(); i++){
        cout << test_1[0].coord[i] << ' ';
    }
    cout << '\n' << "exp: " << test_1[1].exp << '\n';
    cout << "quantum numbers " << '\n';
    cout << test_1[1].A << '\n';

    return 0;
}
