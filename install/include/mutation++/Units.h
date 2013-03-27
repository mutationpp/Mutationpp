#ifndef UTILS_UNITS_H
#define UTILS_UNITS_H

#include <vector>
#include <string>

namespace Mutation {
    namespace Utilities {

class Units
{
public:

    static std::vector<Units> split(const std::string str);

public:

    Units();
    Units(const char []);
    Units(const std::string& str);
    //Units(std::initializer_list<double> exponents);
    Units(double e1, double e2, double e3, double e4, double e5, double e6, double e7); 
    Units(const Units& units);
    
    double convertToBase(const double number);
    double convertTo(const double number, const Units& units);
    double convertTo(const double number, const std::string& units);
    
    Units& operator=(const Units& right);
    
    friend Units operator*(const Units& left, const double right);
    friend Units operator/(const Units& left, const double right);
    friend Units operator*(const Units& left, const Units& right);
    friend Units operator/(const Units& left, const Units& right);
    friend Units operator^(const Units& left, const double power);

private:

    void initializeFromString(const std::string& str);

private:

    double m_factor;
    double m_exponent[7];
};

    } // namespace Utilities
} // namespace Mutation

#endif // UTILS_UNITS_H
