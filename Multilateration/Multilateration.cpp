// Multilateration.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <cassert>



class Point3D
{
private:
    double x, y, z;

public:
    Point3D(double x, double y, double z) : x(x), y(y), z(z) {}

    double getX() {
        return x;
    }

    double getY() {
        return y;
    }

    double getZ() {
        return z;
    }

};



bool quadraticEquation(double a, double b, double c, std::vector<double>& solution)
{
    double d;
    d = b * b - 4 * a * c;
    if (d > 0)
    {
        solution.push_back((-b + sqrt(d)) / (2 * a));
        solution.push_back((-b - sqrt(d)) / (2 * a));
        return true;
    }
    if (d == 0)
    {
        solution.push_back(-b / 2 * a);
        return true;
    }
    return false;
}



// TODO think if d1, d2, d3 should be double but not int
// TODO think about return type, maybe pair<Point3D, Point3D> is better ?

/// Function finding pair of eventual coordinates of the source by given 3points and their distances to the source
/// \param d1 is distance from point p1 to the source
/// \param d2, d3 distances respectively from p2 and d3 to the source
/// \return Pair of 3D Points are returned as a eventual coordinates of the source (4th point) 
std::vector<Point3D> findFourthPoint(Point3D& p1, double d1, Point3D& p2, double d2, Point3D& p3, double d3) {

    std::vector<Point3D> sound_source_coords;     // result coordinates

    double x1 = p1.getX(),
        y1 = p1.getY(),
        z1 = p1.getZ();
    double x2 = p2.getX(),
        y2 = p2.getY(),
        z2 = p2.getZ();
    double x3 = p3.getX(),
        y3 = p3.getY(),
        z3 = p3.getZ();

    double a11 = 2 * (x3 - x2),
        a12 = 2 * (y3 - y2),
        a13 = 2 * (z3 - z2),
        a21 = 2 * (x1 - x2),
        a22 = 2 * (y1 - y2),
        a23 = 2 * (z1 - z2)
        ;

    double b1 =
        pow(d2, 2) - pow(d3, 2) + pow(x3, 2) - pow(x2, 2) + pow(y3, 2) - pow(y2, 2) + pow(z3, 2) - pow(z2, 2),
        b2 =
        pow(d2, 2) - pow(d1, 2) + pow(x1, 2) - pow(x2, 2) + pow(y1, 2) - pow(y2, 2) + pow(z1, 2) - pow(z2, 2),
        b3 =
        pow(d1, 2) - pow(d3, 2) + pow(x3, 2) - pow(x1, 2) + pow(y3, 2) - pow(y1, 2) + pow(z3, 2) - pow(z1, 2);

    /// following assert must be != 0 because of forbidden dividing by 0
    assert(a21 * a12 - a22 * a11 != 0);
    assert(a22 * a11 - a21 * a12 != 0);
    double A = (a12 * b2 - a22 * b1) / (a21 * a12 - a22 * a11);
    double B = (a23 * a12 - a22 * a13) / (a21 * a12 - a22 * a11);
    double C = (b2 * a11 - a21 * b1) / (a22 * a11 - a21 * a12);
    double D = (a23 * a11 - a21 * a13) / (a22 * a11 - a21 * a12);

    // std::cout << pow(2*B(x2-A) + 2*D(y1-C) - 2*z1, 2) - 4 * (pow(B, 2) + pow(D, 2) + 1)*()

    std::vector<double> solution;

    quadraticEquation(
        pow(B, 2) + pow(D, 2) + 1,
        2 * B * (x1 - A) + 2 * D * (y1 - C) - 2 * z1,
        pow(x1 - A, 2) + pow(y1 - C, 2) + pow(z1, 2) - pow(d1, 2),
        solution
    );

    for (auto z_coord : solution) {
        double x_coord = A - z_coord * B;
        double y_coord = C - z_coord * D;

        sound_source_coords.push_back(
            Point3D(
                x_coord,
                y_coord,
                z_coord
            )
        );
    }

    return sound_source_coords;
}



// TODO change parameters (Point3D --> MicrophoneData)

/// Finds source coords (pair of coords)
/// \param d1 is distance from first microphone to the source
/// \return 2 Points in 3D for a particular distance
std::vector<Point3D> solveThreeMicrophonesProblem(Point3D& p1, Point3D& p2, Point3D& p3, double d1) {
    double delay2 = 2;                            // delay(p1, p2)
    double delay3 = 3;                            // delay(p1, p3)

    double d2 = d1 + delay2;
    double d3 = d1 + delay3;

    return findFourthPoint(p1, d1
                         , p2, d2
                         , p3, d3);
}


/// TODO ISOLATE solving non-linear system with 4 unknown variables
/// Solve for three microphones and eventual distance of the sound source
///  and then check with results for x,y,z if they are solution for fourth equation
void solveFourMicrophonesProblem(Point3D& p1, Point3D& p2, Point3D& p3, Point3D& p4) {
    double fourth_equation_result_1;
    double fourth_equation_result_2;
    double x4 = p4.getX(),
        y4 = p4.getY(),
        z4 = p4.getZ();
    double delay4 = 1.5;  /// delay(p1, p4)
    std::vector<Point3D> threeMicSolutionCoords_1 = solveThreeMicrophonesProblem(p1, p2, p3, 0.373);   /// 4 should be modified
    std::vector<Point3D> threeMicSolutionCoords_2;

    // Point3D first_solution;                              /// eventual sound source coordinates
    // Point3D second_solution;                             /// eventual sound source coordinates

    for (double eventual_distance = 0.373; eventual_distance <= 5000; ++eventual_distance) {
        threeMicSolutionCoords_2 = solveThreeMicrophonesProblem(p1, p2, p3, eventual_distance);

        // following two if's to check if line is defined on eventual_distance
        if (threeMicSolutionCoords_2.size() == 0) {
            threeMicSolutionCoords_1 = threeMicSolutionCoords_2;
            continue;
        }
        if (threeMicSolutionCoords_1.size() == 0) {
            threeMicSolutionCoords_1 = threeMicSolutionCoords_2;
            continue;
        }


        // std::cout << eventual_distance << std::endl;

        fourth_equation_result_1 =
            pow(threeMicSolutionCoords_1[0].getX() - x4, 2) + pow(threeMicSolutionCoords_1[0].getY() - y4, 2) + pow(threeMicSolutionCoords_1[0].getZ() - z4, 2) - pow(eventual_distance - 1 + delay4, 2);
        
        fourth_equation_result_2 =
            pow(threeMicSolutionCoords_2[0].getX() - x4, 2) + pow(threeMicSolutionCoords_2[0].getY() - y4, 2) + pow(threeMicSolutionCoords_2[0].getZ() - z4, 2) - pow(eventual_distance + delay4, 2);
        
        if (fourth_equation_result_1 <= 0 && fourth_equation_result_2 >= 0)
            std::cout << eventual_distance << "\n";
        if (fourth_equation_result_1 >= 0 && fourth_equation_result_2 <= 0)
            std::cout << eventual_distance << "\n";
        // std::cout << fourth_equation_result_1 << " " << fourth_equation_result_2 << "\n";


        fourth_equation_result_1 =
            pow(threeMicSolutionCoords_1[1].getX() - x4, 2) + pow(threeMicSolutionCoords_1[1].getY() - y4, 2) + pow(threeMicSolutionCoords_1[1].getZ() - z4, 2) - pow(eventual_distance - 1 + delay4, 2);
        fourth_equation_result_2 =
            pow(threeMicSolutionCoords_2[1].getX() - x4, 2) + pow(threeMicSolutionCoords_2[1].getY() - y4, 2) + pow(threeMicSolutionCoords_2[1].getZ() - z4, 2) - pow(eventual_distance + delay4, 2);
        
        if (fourth_equation_result_1 <= 0 && fourth_equation_result_2 >= 0)
            std::cout << eventual_distance << "\n";
        if (fourth_equation_result_1 >= 0 && fourth_equation_result_2 <= 0)
            std::cout << eventual_distance << "\n";
        // std::cout << fourth_equation_result_1 << " " << fourth_equation_result_2 << "\n";

        threeMicSolutionCoords_1 = threeMicSolutionCoords_2;
    }

}




int main()
{
    double d = 38;
    // microphones coordinates
    Point3D p1(-1, 6, 4),
        p2(3, 6, 3),
        p3(-8, 3, 1),
        p4(-5, 2, 14);
    
    for (auto i : solveThreeMicrophonesProblem(p1, p2, p3, d))
        std::cout << i.getX() << " " << i.getY() << " " << i.getZ() << std::endl;
    
    solveFourMicrophonesProblem(p1, p2, p3, p4);
}
