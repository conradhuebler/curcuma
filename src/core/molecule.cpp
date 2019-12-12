#include "molecule.h"
#include <cstdio>
#include <fstream>
#include <sstream>
#include <istream>
#include <iostream>
#include <cstring>
#include <cmath>
#include <array>



void Molecule::print_geom() const
{
    std::cout << natom  << endl << endl;
    for(int i=0; i < natom; i++)
    {
        if(atoms[i].size())
            printf("%s %8.5f %8.5f %8.5f\n", atoms[i].c_str(), geom[i][0], geom[i][1], geom[i][2]);
    }
}

void Molecule::printAtom(int i) const
{
    --i;
    if(i < natom)
        printf("%s %8.5f %8.5f %8.5f", atoms[i].c_str(), geom[i][0], geom[i][1], geom[i][2]);
}


void Molecule::translate(double x, double y, double z)
{
    for(int i=0; i < natom; i++) {
        geom[i][0] += x;
        geom[i][1] += y;
        geom[i][2] += z;
    }
}

Molecule::Molecule(int n, int q)
{
    natom = n;
    charge = q;
    
    atoms = std::vector<std::string>(natom);
    for(int i=0; i < natom; i++)
    {
        std::array<double, 3> atom = {0,0,0};
        geom.push_back(  atom );
//         geom[i][0] = 0;
//         geom[i][1] = 0;
//         geom[i][2] = 0;
//         zvals[i] = 0;
    }
}
Molecule::~Molecule()
{ 

}


double Molecule::Distance(int i, int j) const
{
    if(i >= natom || j >= natom)
        return 0;
    
    double x_i = geom[i][0];
    double x_j = geom[j][0];
    
    double y_i = geom[i][1];
    double y_j = geom[j][1];
    
    double z_i = geom[i][2];
    double z_j = geom[j][2];
    
    return sqrt( ( ((x_i-x_j)*(x_i-x_j))+((y_i-y_j)*(y_i-y_j))+((z_i-z_j)*(z_i-z_j)) ));
    
}

double Molecule::DotProduct(std::array<double, 3> pos1, std::array<double, 3> pos2) const
{
    return pos1[0]*pos2[0]+pos1[1]*pos2[1]+pos1[2]*pos2[2];
}


double Molecule::angle(int atom1, int atom2, int atom3) const
{

    std::array<double, 3> atom_0 = { geom[atom2 - 1] }; // Proton
    std::array<double, 3> atom_1 = { geom[atom3 - 1] }; // Acceptor
    std::array<double, 3> atom_2 = { geom[atom1 - 1] }; // Donor
    
    
    std::array<double, 3> vec_1 = { atom_0[0]-atom_1[0], atom_0[1]-atom_1[1], atom_0[2]-atom_1[2] };
    std::array<double, 3> vec_2 = { atom_0[0]-atom_2[0], atom_0[1]-atom_2[1], atom_0[2]-atom_2[2] };
    
    return acos(DotProduct(vec_1,vec_2)/(sqrt(DotProduct(vec_1,vec_1)*DotProduct(vec_2,vec_2))))*360/2.0/pi;
}



void Molecule::setAtom(const std::string& internal, int i)
{
    std::vector<std::string > elements;
    std::string element;
    const char *delim = " ";
    for (const char & c : internal)
    {
        if(*delim != c)
            element.push_back(c);
        else
        {
            if(element.size())
                elements.push_back(element);
            element.clear();
        }
    }
                    
    elements.push_back(element);

    if(elements.size() == 7)
    {
        atoms[i] = elements[0];
        
        int atom_1 = stoi(elements[1]);
        int atom_2 = stoi(elements[2]);
        int atom_3 = stoi(elements[3]);
        
        double r_ij = stod(elements[4]);
        double omega = stod(elements[5]);
        double theta = stod(elements[6]);
        
        double r_ik = Distance(atom_2, atom_1);
        
        double x = r_ik + r_ij*cos(omega);
        double y = r_ij*sin(omega);
        double z = 0;
        geom[i][0] = x;
        geom[i][1] = y;
        geom[i][2] = z;
    }
        
}

void Molecule::setXYZ(const std::string& internal, int i)
{
    std::vector<std::string > elements;
    std::string element;
    const char *delim = " ";
    for (const char & c : internal)
    {
        if(*delim != c)
            element.push_back(c);
        else
        {
            if(element.size())
                elements.push_back(element);
            element.clear();
        }
    }

    elements.push_back(element);

    if(elements.size() >= 4)
    {
        atoms[i] = elements[0];
        
        double x = stod(elements[1]);
        double y = stod(elements[2]);
        double z = stod(elements[3]);
        
        geom[i][0] = x;
        geom[i][1] = y;
        geom[i][2] = z;
    }
        
}

Geometry Molecule::getGeometry() const
{
    Geometry geometry(geom.size(), 3);

    for(int i = 0; i < geom.size(); ++i)
    {
        geometry(i, 0) = geom[i][0];
        geometry(i, 1) = geom[i][1];
        geometry(i, 2) = geom[i][2];
    }
    return geometry;
}

bool Molecule::setGeometry(const Geometry &geometry)
{
    if(geometry.rows() != geom.size())
        return false;

    for(int i = 0; i < geom.size(); ++i)
    {
        geom[i][0] = geometry(i, 0);
        geom[i][1] = geometry(i, 1);
        geom[i][2] = geometry(i, 2);
    }

    return true;
}

Position  Molecule::Centroid(bool hydrogen) const
{
    Position position = {0 , 0 , 0};
    Geometry geom = getGeometry();
    if(hydrogen)
    {
        for(int i = 0; i < geom.rows(); ++i)
        {
            position += geom.row(i);
        }

        position /= double(geom.rows());
    }else
    {
        throw -1;
    }
    return position;
}
