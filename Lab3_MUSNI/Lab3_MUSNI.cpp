#include <iostream>
#include <cmath>
using namespace std;

void simplify(double matrix[10][11], int var);



int main() {
    int var, choice;
    double matrix[10][11], solution[10];;
    char opt;
    
    do {
        menu:
        cout << "DIRECT METHOD FOR SOLVING LINEAR SYSTEM\n\n"
            << "\t1.] Gauss Elimination Method\n"
            << "\t2.] Gauss Elimination with Maximum Pivot Strategy\n"
            << "\t3.] Gauss Jordan Method\n\n"
            << "Enter choice: ";
        cin >> choice;
        
        switch (choice) {
        case 1:
            system("cls");
            one:
            cout << "Gauss Elimination Method\n\n"
                << "Enter the number of variables (up to 10): ";
            cin >> var;

            if (cin.fail() || var < 1 || var > 10) {
                cout << "Invalid input! Please enter an integer between 1 and 10.\n";
                cin.clear();
                cin.ignore(1000, '\n');
                system("pause");
                system("cls");
                goto one;
            }

            cout << "\nEnter the augmented matrix:\n";
            for (int i = 0; i < var; i++) {
                for (int j = 0; j <= var; j++) {
                    cout << "[" << i << "][" << j << "]: ";
                    cin >> matrix[i][j];
                }
            }
            
            for (int i = 0; i < var; i++) {
                int maxRow = i;
                for (int k = i + 1; k < var; k++) {
                    if (fabs(matrix[k][i]) > fabs(matrix[maxRow][i])) {
                        maxRow = k;
                    }
                }
                
                if (maxRow != i) {
                    for (int j = 0; j <= var; j++) {
                        swap(matrix[i][j], matrix[maxRow][j]);
                    }
                }
                
                double pivot = matrix[i][i];
                if (pivot == 0) {
                    cout << "Mathematical Error: Zero pivot encountered.\n";
                    return 0;
                }
                
                for (int j = i; j <= var; j++) {
                    matrix[i][j] /= pivot;
                }
                
                for (int k = i + 1; k < var; k++) {
                    double factor = matrix[k][i];
                    for (int j = i; j <= var; j++) {
                        matrix[k][j] -= factor * matrix[i][j];
                    }
                }
            }
            
            for (int i = var - 1; i >= 0; i--) {
                double sum = matrix[i][var];
                for (int j = i + 1; j < var; j++) {
                    sum -= matrix[i][j] * solution[j];
                }
                solution[i] = sum; 
            }
            
            cout << "\nFinal Upper Triangular Matrix:\n";
            for (int i = 0; i < var; i++) {
                for (int j = 0; j <= var; j++) {
                    cout << matrix[i][j] << "\t\t";
                }
                cout << "\n\n";
            }

            simplify(matrix, var);
            break;
            
        case 2:
            system("cls");
            two:
            cout << "Gaussian Elimination with Maximum Pivot Strategy\n\n"
                << "Enter the number of variables (up to 10): ";
            cin >> var;

            if (cin.fail() || var < 1 || var > 10) {
                cout << "Invalid input! Please enter an integer between 1 and 10.\n";
                cin.clear();
                cin.ignore(1000, '\n');
                system("pause");
                system("cls");
                goto two;
            }

            cout << "\nEnter the augmented matrix:\n";
            for (int i = 0; i < var; i++) {
                for (int j = 0; j <= var; j++) {
                    cout << "[" << i << "][" << j << "]: ";
                    cin >> matrix[i][j];
                }
            }
            
            for (int i = 0; i < var - 1; i++) {
                int maxRow = i;
                for (int k = i + 1; k < var; k++) {
                    if (fabs(matrix[k][i]) > fabs(matrix[maxRow][i])) {
                        maxRow = k;
                    }
                }
                
                if (maxRow != i) {
                    for (int j = 0; j <= var; j++) {
                        swap(matrix[i][j], matrix[maxRow][j]);
                    }
                }
                
                if (matrix[i][i] == 0) {
                    cout << "Mathematical Error: Zero pivot encountered.\n";
                    return 0;
                }
                
                for (int k = i + 1; k < var; k++) {
                    double factor = matrix[k][i] / matrix[i][i];
                    for (int j = i; j <= var; j++) {
                        matrix[k][j] -= factor * matrix[i][j];
                    }
                }
            }
            
            for (int i = var - 1; i >= 0; i--) {
                double sum = matrix[i][var];
                for (int j = i + 1; j < var; j++) {
                    sum -= matrix[i][j] * solution[j];
                }
                solution[i] = sum / matrix[i][i];
            }
            
            cout << "\nFinal Upper Triangular Matrix:\n";
            for (int i = 0; i < var; i++) {
                for (int j = 0; j <= var; j++) {
                    cout << matrix[i][j] << "\t\t";
                }
                cout << "\n\n";
            }
            
            simplify(matrix, var);
            break;
            
        case 3:
            system("cls");
            three:
            cout << "Gauss Jordan Method of Solving Linear Equations\n\n"
                << "Enter the number of variables (up to 10): ";
            cin >> var;

            if (cin.fail() || var < 1 || var > 10) {
                cout << "Invalid input! Please enter an integer between 1 and 10.\n";
                cin.clear();
                cin.ignore(1000, '\n');
                system("pause");
                system("cls");
                goto three;
            }
            
            cout << "\nEnter the augmented matrix:\n";
            for (int i = 0; i < var; i++) {
                for (int j = 0; j <= var; j++) {
                    cout << "[" << i << "][" << j << "]: ";
                    cin >> matrix[i][j];
                }
            }

            for (int i = 0; i < var; i++) {
                int pivotRow = i;
                for (int k = i + 1; k < var; k++) {
                    if (fabs(matrix[k][i]) > fabs(matrix[pivotRow][i])) {
                        pivotRow = k;
                    }
                }

                if (pivotRow != i) {
                    for (int j = 0; j <= var; j++) {
                        swap(matrix[i][j], matrix[pivotRow][j]);
                    }
                }

                double pivot = matrix[i][i];
                if (pivot == 0) {
                    cout << "Mathematical Error: Zero pivot encountered.\n";
                    return 0;
                }
                for (int j = 0; j <= var; j++) {
                    matrix[i][j] /= pivot;
                }

                for (int k = 0; k < var; k++) {
                    if (k != i) {
                        double factor = matrix[k][i];
                        for (int j = 0; j <= var; j++) {
                            matrix[k][j] -= factor * matrix[i][j];
                        }
                    }
                }
            }

            cout << "\nFinal RREF Matrix:\n";
            for (int i = 0; i < var; i++) {
                for (int j = 0; j <= var; j++) {
                    cout << matrix[i][j] << "\t\t";
                }
                cout << "\n\n";
            }
            break;
            
        default:
            cout << "Invalid choice! Please choose from 1 to 3 only.\n";
            cin.clear();
            cin.ignore(1000, '\n');
            system("pause");
            system("cls");
            goto menu;
            break;
        }
        
        cout << "\nSolution:\n";
        for (int i = 0; i < var; i++) {
            cout << "x" << i + 1 << " = " << matrix[i][var] << endl;
        }
           
        cout << "\nWould you like to try again? (Press Y/y to continue...) ";
        cin >> opt;
        system("cls");
        
    } while ( opt == 'y' || opt == 'Y');
    

    return 0;
}










void simplify(double matrix[10][11], int var) {
    for (int i = 0; i < var; i++) {
        int pivotRow = i;

        for (int k = i + 1; k < var; k++) {
            if (fabs(matrix[k][i]) > fabs(matrix[pivotRow][i])) {
                pivotRow = k;
            }
        }

        if (pivotRow != i) {
            for (int j = 0; j <= var; j++) {
                swap(matrix[i][j], matrix[pivotRow][j]);
            }
        }

        double pivot = matrix[i][i];
        if (pivot == 0) {
            cout << "Mathematical Error: Zero pivot encountered.\n";
            return;
        }

        for (int j = 0; j <= var; j++) {
            matrix[i][j] /= pivot;
        }

        for (int k = 0; k < var; k++) {
            if (k != i) {
                double factor = matrix[k][i];
                for (int j = 0; j <= var; j++) {
                    matrix[k][j] -= factor * matrix[i][j];
                }
            }
        }
    }
}
