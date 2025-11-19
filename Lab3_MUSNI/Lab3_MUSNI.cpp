#include <iostream>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <limits>
using namespace std;

void printMatrix(double a[10][11], int n) {
    cout << fixed << setprecision(6);
    for (int i = 0; i < n; i++) {
        cout << "| ";
        for (int j = 0; j <= n; j++) {
            cout << setw(12) << a[i][j] << " ";
        }
        cout << "|\n";
    }
    cout << endl;
}

void displayEquations(double a[10][11], int n) {
    cout << "\nSystem of Equations:\n";
    for (int i = 0; i < n; i++) {
        cout << "Eq" << i + 1 << ":  ";
        for (int j = 0; j < n; j++) {
            if(a[i][j] == 1)
                cout << "x" << j + 1;
            else
                cout << a[i][j] << "x" << j + 1;
            if (j < n - 1) cout << "  +  ";
        }
        cout << "  =  " << a[i][n] << "\n";
    }
    cout << endl;
}

void displayIterationFormulas(double a[10][11], int n) {
    cout << "\nIteration Formulas (for Seidel / Jacobi):\n\n";

    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = ( " << a[i][n];

        // subtract terms
        for (int j = 0; j < n; j++) {
            if (j != i) {
                if (a[i][j] >= 0) cout << " - " << a[i][j] << "x" << j + 1;
                else cout << " + " << fabs(a[i][j]) << "x" << j + 1;
            }
        }

        cout << " ) / " << a[i][i];
        cout << "\n\n";
    }
}

void swapRows(double a[10][11], int r1, int r2, int n) {
    for (int j = 0; j <= n; j++) {
        double temp = a[r1][j];
        a[r1][j] = a[r2][j];
        a[r2][j] = temp;
    }
}

void gaussElimination(double a[10][11], int n) {
    double sol[10];
    
    displayEquations(a, n);

    cout << "\nInitial Matrix:\n";
    printMatrix(a, n);

    for (int k = 0; k < n - 1; k++) {

        if (fabs(a[k][k]) < 1e-12) {
            int sw = -1;
            for (int r = k + 1; r < n; r++) {
                if (fabs(a[r][k]) > 1e-12) {
                    sw = r;
                    break;
                }
            }
            if (sw != -1) {
                cout << "R" << k + 1 << " <-> R" << sw + 1 << "\n";
                swapRows(a, k, sw, n);
                printMatrix(a, n);
            }
            else {
                cout << "Zero pivot encountered.\n";
                return;
            }
        }

        for (int i = k + 1; i < n; i++) {
            double factor = a[i][k] / a[k][k];
            cout << "R" << i + 1 << " = R" << i + 1 << " - (" << factor << ")R" << k + 1 << "\n";

            for (int j = k; j <= n; j++) {
                a[i][j] -= factor * a[k][j];
            }
            printMatrix(a, n);
        }
    }

    cout << "Normalize diagonal rows:\n";
    for (int i = 0; i < n; i++) {
        double pivot = a[i][i];
        cout << "R" << i + 1 << " = R" << i + 1 << " / (" << pivot << ")\n";
        for (int j = i; j <= n; j++) {
            a[i][j] /= pivot;
        }
        printMatrix(a, n);
    }

    cout << "Back substitution:\n";
    for (int i = n - 1; i >= 0; i--) {
        double s = a[i][n];
        cout << "x" << i + 1 << " = " << s;
        for (int j = i + 1; j < n; j++) {
            cout << " - (" << a[i][j] << ")*x" << j + 1;
            s -= a[i][j] * sol[j];
        }
        sol[i] = s;
        cout << "  = " << sol[i] << "\n";
    }

    cout << "\nSolution:\n";
    for (int i = 0; i < n; i++) cout << "x" << i + 1 << " = " << sol[i] << endl;
}

void gemps(double a[10][11], int n) {
    displayEquations(a, n);

    cout << "\nInitial Matrix:\n";
    printMatrix(a, n);

    for (int k = 0; k < n; k++) {

        int maxRow = k;
        for (int r = k + 1; r < n; r++) {
            if (fabs(a[r][k]) > fabs(a[maxRow][k])) {
                maxRow = r;
            }
        }
        if (maxRow != k) {
            cout << "R" << k + 1 << " <-> R" << maxRow + 1 << "\n";
            swapRows(a, k, maxRow, n);
            printMatrix(a, n);
        }

        double pivot = a[k][k];
        cout << "R" << k + 1 << " = R" << k + 1 << " / (" << pivot << ")\n";
        for (int j = k; j <= n; j++) {
            a[k][j] /= pivot;
        }
        printMatrix(a, n);

        for (int i = k + 1; i < n; i++) {
            double factor = a[i][k];
            cout << "R" << i + 1 << " = R" << i + 1 << " - (" << factor << ")R" << k + 1 << "\n";
            for (int j = k; j <= n; j++) {
                a[i][j] -= factor * a[k][j];
            }
            printMatrix(a, n);
        }
    }

    cout << "Backward substitution:\n";
    double sol[10];
    for (int i = n - 1; i >= 0; i--) {
        double val = a[i][n];
        cout << "x" << i + 1 << " = " << val;
        for (int j = i + 1; j < n; j++) {
            cout << " - (" << a[i][j] << ")*x" << j + 1;
            val -= a[i][j] * sol[j];
        }
        sol[i] = val;
        cout << " = " << sol[i] << "\n";
    }

    cout << "\nSolution:\n";
    for (int i = 0; i < n; i++) cout << "x" << i + 1 << " = " << sol[i] << endl;
}

void gaussJordan(double a[10][11], int n) {
    displayEquations(a, n);

    cout << "\nInitial Matrix:\n";
    printMatrix(a, n);

    for (int k = 0; k < n; k++) {

        int pivotRow = k;
        for (int r = k + 1; r < n; r++) {
            if (fabs(a[r][k]) > fabs(a[pivotRow][k])) {
                pivotRow = r;
            }
        }
        if (pivotRow != k) {
            cout << "R" << k + 1 << " <-> R" << pivotRow + 1 << "\n";
            swapRows(a, k, pivotRow, n);
            printMatrix(a, n);
        }

        double pivot = a[k][k];
        cout << "R" << k + 1 << " = R" << k + 1 << " / (" << pivot << ")\n";
        for (int j = 0; j <= n; j++) {
            a[k][j] /= pivot;
        }
        printMatrix(a, n);

        for (int i = 0; i < n; i++) {
            if (i != k) {
                double factor = a[i][k];
                cout << "R" << i + 1 << " = R" << i + 1 << " - (" << factor << ")R" << k + 1 << "\n";
                for (int j = 0; j <= n; j++) {
                    a[i][j] -= factor * a[k][j];
                }
                printMatrix(a, n);
            }
        }
    }

    cout << "Solution (RREF):\n";
    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = " << a[i][n] << endl;
    }
}

void gaussSeidel(double a[10][11], int n) {

    double x[10], prev_x[10], e[10];
    
    displayEquations(a, n);
    displayIterationFormulas(a, n);

    for (int i = 0; i < n; i++) {
        x[i] = 0;
        prev_x[i] = 0;
    }

    cout << fixed << setprecision(6);
    cout << "\n--------------------------------------------------------------------------------------------------------------------------------------\n";
    cout << "k\t";
    for (int i = 0; i < n; i++) cout << "x" << i + 1 << "\t\t\t";
    for (int i = 0; i < n; i++) cout << "|e(x" << i + 1 << ")|\t\t\t";
    cout << "\n--------------------------------------------------------------------------------------------------------------------------------------\n";

    cout << "0\t";
    for (int i = 0; i < n; i++) cout << x[i] << "\t\t";
    for (int i = 0; i < n; i++) cout << "-\t\t\t";
    cout << "\n";

    const double TOL = 0.001;

    bool solved = false;

    for (int k = 1; k <= 15; k++) {

        for (int i = 0; i < n; i++) prev_x[i] = x[i];

        for (int i = 0; i < n; i++) {
            double sum = a[i][n];
            for (int j = 0; j < n; j++) {
                if (j != i)
                    sum -= a[i][j] * x[j];
            }
            x[i] = sum / a[i][i];
        }

        bool allSatisfied = true;

        for (int i = 0; i < n; i++) {
            e[i] = fabs(x[i] - prev_x[i]);
            if (e[i] >= TOL) allSatisfied = false;
        }

        cout << k << "\t";
        for (int i = 0; i < n; i++) cout << x[i] << "\t\t";
        for (int i = 0; i < n; i++) cout << e[i] << "\t\t";
        cout << "\n";

        if (allSatisfied) {
            cout << "\nAll |e(k)| < 0.001 satisfied at iteration " << k << "\n";
            solved = true;
            break;
        }
    }

    if (!solved) {
        cout << "\nGauss-Seidel cannot solve the system (failed to converge after 15 iterations).\n";
        return;
    }

    cout << "\nFinal Approximation:\n";
    for (int i = 0; i < n; i++) cout << "x" << i + 1 << " = " << x[i] << endl;
}

void jacobi(double a[10][11], int n) {
    double x[10], x_new[10], e[10];
    
    displayEquations(a, n);
    displayIterationFormulas(a, n);

    for (int i = 0; i < n; i++) {
        x[i] = 0;
        x_new[i] = 0;
    }

    cout << fixed << setprecision(6);
    cout << "\n--------------------------------------------------------------------------------------------------------------------------------------\n";
    cout << "k\t";
    for (int i = 0; i < n; i++) cout << "x" << i + 1 << "\t\t\t";
    for (int i = 0; i < n; i++) cout << "|e(x" << i + 1 << ")|\t\t\t";
    cout << "\n--------------------------------------------------------------------------------------------------------------------------------------\n";

    cout << "0\t";
    for (int i = 0; i < n; i++) cout << x[i] << "\t\t";
    for (int i = 0; i < n; i++) cout << "-\t\t\t";
    cout << "\n";

    const double TOL = 0.001;

    bool solved = false;

    for (int k = 1; k <= 15; k++) {

        for (int i = 0; i < n; i++) {
            double sum = a[i][n];
            for (int j = 0; j < n; j++) {
                if (j != i)
                    sum -= a[i][j] * x[j];
            }
            x_new[i] = sum / a[i][i];
        }

        bool allSatisfied = true;

        for (int i = 0; i < n; i++) {
            e[i] = fabs(x_new[i] - x[i]);
            if (e[i] >= TOL) allSatisfied = false;
        }

        cout << k << "\t";
        for (int i = 0; i < n; i++) cout << x_new[i] << "\t\t";
        for (int i = 0; i < n; i++) cout << e[i] << "\t\t";
        cout << "\n";

        if (allSatisfied) {
            cout << "\nAll |e(k)| < 0.001 satisfied at iteration " << k << "\n";
            solved = true;
            break;
        }

        for (int i = 0; i < n; i++) x[i] = x_new[i];
    }

    if (!solved) {
        cout << "\nJacobi cannot solve the system (failed to converge after 15 iterations).\n";
        return;
    }

    cout << "\nFinal Approximation:\n";
    for (int i = 0; i < n; i++) cout << "x" << i + 1 << " = " << x[i] << endl;
}

int main() {
    int n, choice;
    double a[10][11];
    char opt;

    do {
    menu:
        cout << "DIRECT & ITERATIVE METHODS FOR SOLVING LINEAR SYSTEM\n\n";
        cout << "\t1.] Gauss Elimination Method\n";
        cout << "\t2.] Gauss Elimination with Maximum Pivot Strategy\n";
        cout << "\t3.] Gauss Jordan Method\n";
        cout << "\t4.] Gauss Seidel Method\n";
        cout << "\t5.] Gauss Jacobi Method\n\n";
        cout << "Enter choice: ";
        cin >> choice;

        switch (choice) {
        case 1:
            system("cls");
        one:
            cout << "Gauss Elimination Method\n\n"
                << "Enter the number of variables (up to 10): ";
            cin >> n;

            if (cin.fail() || n < 1 || n > 10) {
                cout << "Invalid input! Please enter an integer between 1 and 10.\n";
                cin.clear();
                cin.ignore(1000, '\n');
                system("pause");
                system("cls");
                goto one;
            }

            cout << "\nEnter the augmented matrix:\n";
            for (int i = 0; i < n; i++) {
                for (int j = 0; j <= n; j++) {
                    string input;
                    double value;
                    while (true) {
                        cout << "[" << i << "][" << j << "]: ";
                        if (i == 0 && j == 0) {
                            cin >> input;
                            cin.ignore(1000, '\n');
                        }
                        else {
                            getline(cin, input);
                        }
                        
                        stringstream ss(input);
                        if (ss >> value && ss.eof()) { 
                            a[i][j] = value;
                            break;
                        }
                        else {
                            cout << "Invalid input! Please enter a number.\n";
                        }
                    }
                }
            }

            gaussElimination(a, n);
            break;
        case 2:
            system("cls");
        two:
            cout << "Gaussian Elimination with Maximum Pivot Strategy\n\n"
                << "Enter the number of variables (up to 10): ";
            cin >> n;

            if (cin.fail() || n < 1 || n > 10) {
                cout << "Invalid input! Please enter an integer between 1 and 10.\n";
                cin.clear();
                cin.ignore(1000, '\n');
                system("pause");
                system("cls");
                goto two;
            }

            cout << "\nEnter the augmented matrix:\n";
            for (int i = 0; i < n; i++) {
                for (int j = 0; j <= n; j++) {
                    string input;
                    double value;
                    while (true) {
                        cout << "[" << i << "][" << j << "]: ";
                        if (i == 0 && j == 0) {
                            cin >> input;
                            cin.ignore(1000, '\n');
                        }
                        else {
                            getline(cin, input);
                        }

                        stringstream ss(input);
                        if (ss >> value && ss.eof()) { 
                            a[i][j] = value;
                            break; 
                        }
                        else {
                            cout << "Invalid input! Please enter a number.\n";
                        }
                    }
                }
            }

            gemps(a, n);
            break;
        case 3:
            system("cls");
        three:
            cout << "Gauss Jordan Method\n\n"
                << "Enter the number of variables (up to 10): ";
            cin >> n;

            if (cin.fail() || n < 1 || n > 10) {
                cout << "Invalid input! Please enter an integer between 1 and 10.\n";
                cin.clear();
                cin.ignore(1000, '\n');
                system("pause");
                system("cls");
                goto three;
            }

            cout << "\nEnter the augmented matrix:\n";
            for (int i = 0; i < n; i++) {
                for (int j = 0; j <= n; j++) {
                    string input;
                    double value;
                    while (true) {
                        cout << "[" << i << "][" << j << "]: ";
                        if (i == 0 && j == 0) {
                            cin >> input;
                            cin.ignore(1000, '\n');
                        }
                        else {
                            getline(cin, input);
                        }

                        stringstream ss(input);
                        if (ss >> value && ss.eof()) {
                            a[i][j] = value;
                            break;
                        }
                        else {
                            cout << "Invalid input! Please enter a number.\n";
                        }
                    }
                }
            }

            gaussJordan(a, n);
            break;
        case 4:
            system("cls");
        four:
            cout << "Gauss Seidel Method\n\n"
                << "Enter the number of variables (up to 10): ";
            cin >> n;

            if (cin.fail() || n < 1 || n > 10) {
                cout << "Invalid input! Please enter an integer between 1 and 10.\n";
                cin.clear();
                cin.ignore(1000, '\n');
                system("pause");
                system("cls");
                goto four;
            }

            cout << "\nEnter the augmented matrix:\n";
            for (int i = 0; i < n; i++) {
                for (int j = 0; j <= n; j++) {
                    string input;
                    double value;
                    while (true) {
                        cout << "[" << i << "][" << j << "]: ";
                        if (i == 0 && j == 0) {
                            cin >> input;
                            cin.ignore(1000, '\n');
                        }
                        else {
                            getline(cin, input);
                        }

                        stringstream ss(input);
                        if (ss >> value && ss.eof()) {
                            a[i][j] = value;
                            break;
                        }
                        else {
                            cout << "Invalid input! Please enter a number.\n";
                        }
                    }
                }
            }
            
            gaussSeidel(a, n);
            break;
        case 5:
            system("cls");
        five:
            cout << "Gauss Jacobi Method\n\n"
                << "Enter the number of variables (up to 10): ";
            cin >> n;

            if (cin.fail() || n < 1 || n > 10) {
                cout << "Invalid input! Please enter an integer between 1 and 10.\n";
                cin.clear();
                cin.ignore(1000, '\n');
                system("pause");
                system("cls");
                goto five;
            }

            cout << "\nEnter the augmented matrix:\n";
            for (int i = 0; i < n; i++) {
                for (int j = 0; j <= n; j++) {
                    string input;
                    double value;
                    while (true) {
                        cout << "[" << i << "][" << j << "]: ";
                        if (i == 0 && j == 0) {
                            cin >> input;
                            cin.ignore(1000, '\n');
                        }
                        else {
                            getline(cin, input);
                        }

                        stringstream ss(input);
                        if (ss >> value && ss.eof()) {
                            a[i][j] = value;
                            break;
                        }
                        else {
                            cout << "Invalid input! Please enter a number.\n";
                        }
                    }
                }
            }
            
            jacobi(a, n);
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

        cout << "\n\nWould you like to try again? (Press Y/y to continue...) ";
        cin >> opt;
        system("cls");

    } while (opt == 'Y' || opt == 'y');

    return 0;
}
