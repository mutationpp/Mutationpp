/**
 * @file LookupTable.h
 * 
 * @brief Defines the InterpolationScheme enum and LookupTable class.
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef LOOKUPTABLE_H
#define LOOKUPTABLE_H

#include <algorithm>
#include <cstring>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <cmath>

#include "Errors.h"
#include "Functors.h"

namespace Mutation {
    namespace Utilities {

/**
 * Enumerates possible interpolation schemes that can be used in the
 * LookupTable::lookup() function.
 *
 * @see LookupTable::lookup()
 */
enum InterpolationScheme {
    NEAREST_INDEX, ///< returns the nearest function values in the table
    LINEAR,        ///< uses linear interpolation scheme
    EXPONENTIAL    ///< uses exponential interpolation
};

/**
 * An efficient lookup table abstract base class.
 *
 * A lookup table consists of a set of one or more functions tabulated against
 * a set of indices (independent value).  A lookup table can be used to tabulate
 * a function when it is more efficient to lookup or interpolate the function
 * values from precomputed values then it is to compute them on the fly.
 *
 * <b>Example usage:</b>
 * @code
 * // Load lookup table from file and create array to hold function values
 * LookupTable<double> table("my_lookup_table.dat");
 * double *functions = new double [table.nFunctions()];
 * 
 * // Lookup the functions' values at x using linear interpolation
 * table.lookup(x, functions, LINEAR);
 * @endcode
 *
 * @note While this class is templated so that it can be used for a variety of
 * index and data types, it was designed and tested with only numerical data 
 * types in mind.  Therefore, care should be taken if this class is to be used 
 * with non-numeric types.
 *
 * @todo Need to add operator= and copy constructor functions to protect the
 * pointer.
 *
 * @see LookupTable::lookup()
 */
template<typename IndexType, typename DataType, typename FunctionType> 
class LookupTable
{
public:

    /**
     * Constructs a LookupTable using data loaded from a file.
     *
     * @param file_name  the name of the file to load.
     */
    explicit LookupTable(const std::string &file_name);
    
    /**
     * Constructs a LookupTable using a prescribed function and error tolerance.
     *
     * The table is created such that the maximum error obtained on any given
     * lookup is below a specified tolerance associated with a specified 
     * interpolation scheme.
     *
     * @param low        the low index for the table range
     * @param high       the high index for the table range
     * @param max_error  the allowed error tolerance associated with scheme
     *                   where 
     *                   \f[ error = \max_i |\frac{interp_i}{exact_i} - 1| \f]
     * @param scheme     the interpolation method to be used to determine the
     *                   maximum error in the table
     */
    LookupTable(
        IndexType low, IndexType high, int nfuncs, 
        typename FunctionType::DataProvider& provider, 
        double max_error = 0.01, InterpolationScheme scheme = LINEAR);
    
    /**
     * Constructs a LookupTable using a prescribed number of rows.
     *
     * Given the number of rows to include in the table and the index range,
     * the table is generated using constant spacing between rows.
     *
     * @param low     the low index for the table range
     * @param high    the high index for the table range
     * @param nrows   the number of rows to include in the table
     */
    LookupTable(
        IndexType low, IndexType high, int nrows, int nfuncs,
        typename FunctionType::DataProvider& provider);
    
    /**
     * Destructor.
     */
    virtual ~LookupTable();
    
    /**
     * Interpolates table data between two table entries nearest a search index.
     *
     * The two nearest data points to the index value in the table are first 
     * found using one of two methods.  If the change in index values between
     * each neighboring row is constant for the entire table, then a simple hash
     * function is used to compute the rows nearest the lookup index directly.
     * If however the change in index values between data entries varies, a 
     * binary search algorithm is employed to find the nearest rows in O(log(n))
     * time where n is the number of rows in the table.
     *
     * Once the two nearest rows are determined, an interpolation scheme is 
     * employeed to compute the table values at the lookup index.  Currently,
     * there are three schemes implemented.
     *
     * - Nearest Index: Chooses the table entry closest to the lookup index.
     * - Linear:
     *     \f[ y = y_1 + \frac{x - x_1}{x_2 - x_1} (y_2 - y_1) \f]
     * - Exponential:
     *     \f[ 
     *         y = \exp(\ln y_1 + \frac{x - x_1}{x_2 - x_1} (\ln y_2 - \ln y_1))
     *     \f]
     *
     * @param index     the independent value index.
     * @param p_values  pointer to T array that will hold function values on
     *                  return.
     * @param scheme    the interpolation scheme to use between table values.
     * 
     * @see InterpolationScheme
     */
    void lookup(
        const IndexType &index, DataType* const p_values, 
        const InterpolationScheme scheme = LINEAR) const
    {
        lookup(
            index, 0, m_num_functions, p_values,
            Mutation::Numerics::Equals<DataType>(), scheme);
    }

    void lookup(
        const IndexType& index, const int start, const int end,
        DataType* const p_values, const InterpolationScheme scheme = LINEAR) const
    {
        lookup(
            index, start, end, p_values,
            Mutation::Numerics::Equals<DataType>(), scheme);
    }

    template <typename OP>
    void lookup(
        const IndexType& index, DataType* const p_values, const OP& op,
        const InterpolationScheme scheme = LINEAR) const
    {
        lookup(
            index, 0, m_num_functions, p_values, op, scheme);
    }

    template <typename OP>
    void lookup(
        const IndexType& index, const int start, const int end,
        DataType* const p_values, const OP& op, const InterpolationScheme scheme)
        const;
         
    
    /**
     * Saves the lookup table to a file.
     *
     * @param file_name  the name of the file to save to, it will be overwritten
     */
    void save(const std::string &file_name) const;
    
    /**
     * Returns number of indices (rows) in the table.
     */
    unsigned int nIndices() const { return m_num_indices; }
    
    /**
     * Returns number of functions (columns) in the table.
     */
    unsigned int nFunctions() const { return m_num_functions; }
    
    /**
     * Returns the first (minimum) index value in the table.
     */
    IndexType minIndex() const { return mp_indices[0]; }
    
    /**
     * Returns the last (maximum) index value in the table.
     */
    IndexType maxIndex() const { return mp_indices[m_num_indices-1]; }

    /**
     * Computes the exact value of each function in this table using the given
     * index value.
     */
    //virtual void tableFunction(IndexType, DataType*) const = 0;

protected:
    
    /**
     * Returns a reference to the data provider for this table.
     */
    typename FunctionType::DataProvider& provider() const { 
        return m_data_provider; 
    }

private:
    
    /**
     * Performs variable intialization common to all constructors.
     */
    void initialize(int nrows, int nfuncs);
    
    /**
     * Loads a lookup table from an ASCII data file.
     *
     * The format of the lookup table file is as follows.
     * 
     * rows functions         \n
     * x1 y1(x1) y2(x1) ...   \n 
     * x2 y1(x2) y2(x2) ...   \n
     *      ...
     *
     * @param file_name  the name of the file to load.
     *
     * @see LookupTable::LookupTable()
     */
    void loadTable(const std::string &file_name);
    
    /**
     * Recurisvely populates a data table represented by a list of 
     * (index, values) pairs until the maximum error in the table is reduced
     * below a given error tolerance using a specified interpolation method.
     *
     * @param table      the list of (index, values) pairs
     * @param high       a list iterator pointing to the end of the desired 
     *                   table range to populate (note only the range between 
     *                   high and --high will be populated)
     * @param func       the function that the lookup table represents
     * @param nfuncs     the number of functions func computes
     * @param max_error  the error tolerance for the table
     * @param scheme     the interpolation method to use to determine error
     *
     * @see LookupTable::LookupTable()
     */
    void populateTable(
        std::list< std::pair<IndexType, DataType *> > &table,
        typename std::list< std::pair<IndexType, DataType *> >::iterator &high,
        int nfuncs, double max_error, InterpolationScheme scheme) const;

private:

    typename FunctionType::DataProvider& m_data_provider;

    unsigned int m_num_indices;
    unsigned int m_num_functions;
    
    size_t m_row_size;
    
    bool m_is_constant_delta;
    
    IndexType *mp_indices;
    DataType  *mp_data;
    
}; // class LookupTable

//==============================================================================

template<typename DataType, typename OP>
inline static void _interpolateNearestIndex(
    const double &ratio, const DataType *const p_y1, const DataType *const p_y2, 
    DataType *const p_out, const int &start, const int &end, const OP& op)
{
    const DataType* const p_values = (ratio < 0.5 ? p_y1 : p_y2);
    for (int i = start; i < end; ++i)
        op(p_out[i], p_values[i]);
}

//==============================================================================

template<typename DataType, typename OP>
inline static void _interpolateLinear(
    const double &ratio, const DataType *const p_y1, const DataType *const p_y2, 
    DataType *const p_out, const int &start, const int &end, const OP& op)
{
    for (int i = start; i < end; ++i)
        op(p_out[i], p_y1[i] + ratio * (p_y2[i] - p_y1[i]));
}

//==============================================================================

template<typename DataType, typename OP>
inline static void _interpolateExponential(
    const double &ratio, const DataType *const p_y1, const DataType *const p_y2, 
    DataType *const p_out, const int &start, const int &end, const OP& op)
{
    DataType log1;
    
    for (int i = start; i < end; ++i) {
        log1 = std::log(p_y1[i]);
        op(p_out[i], exp(log1 + ratio * (std::log(p_y2[i]) - log1)));
    }
}

//==============================================================================

template<typename IndexType, typename DataType, typename FunctionType>
LookupTable<IndexType, DataType, FunctionType>::LookupTable(
    const std::string &file_name)
{
    loadTable(file_name);
}

//==============================================================================

template<typename IndexType, typename DataType, typename FunctionType>
LookupTable<IndexType, DataType, FunctionType>::LookupTable(
    IndexType low, IndexType high, int nfuncs, 
    typename FunctionType::DataProvider& provider, double max_error, 
    InterpolationScheme scheme)
    : m_data_provider(provider)
{
    std::list< std::pair<IndexType, DataType *> > points;
    
    // Initialize the list with the bounding points    
    DataType *p_first = new DataType [nfuncs];
    DataType *p_last  = new DataType [nfuncs];
    
    FunctionType function;
    
    //tableFunction(low, p_first);
    function(low, p_first, provider);
    points.push_back(std::pair<IndexType, DataType *>(low, p_first));
    
    //tableFunction(high, p_last);
    function(high, p_last, provider);
    points.push_back(std::pair<IndexType, DataType *>(high, p_last));
    
    // Now add points to the table in order to recursively reduce error
    populateTable(
        points, --points.end(), nfuncs, max_error, scheme);
    
    // Now we finally have all the points as a list, copy the data into the
    // class members
    initialize(points.size(), nfuncs);
    
    typename std::list< std::pair<IndexType, DataType*> >::iterator iter = 
        points.begin();
    
    for (int index = 0; iter != points.end(); ++iter, ++index) {
        // Copy index
        mp_indices[index] = iter->first;
        
        // Copy function values
        memcpy(mp_data + index * m_num_functions, iter->second, m_row_size);
        
        // Be sure to free the list storge as we go
        delete [] iter->second;
        iter->second = NULL;
    }
    
}

//==============================================================================

template<typename IndexType, typename DataType, typename FunctionType>
LookupTable<IndexType, DataType, FunctionType>::LookupTable(
    IndexType low, IndexType high, int nrows, int nfuncs, 
    typename FunctionType::DataProvider& provider)
    : m_data_provider(provider)
{
    // Initialize table memory
    initialize(nrows, nfuncs);
    m_is_constant_delta = true;
    
    // Compute the index spacing between table rows
    const IndexType delta = (high - low) / static_cast<IndexType>(nrows - 1);
    
    // Fill the table
    FunctionType function;
    
    for (int i = 0; i < m_num_indices; ++i) {
        mp_indices[i] = low + delta * static_cast<IndexType>(i);
        //tableFunction(mp_indices[i], mp_data + i * m_num_functions);
        function(mp_indices[i], mp_data + i * m_num_functions, provider);
    }
}

//==============================================================================

template<typename IndexType, typename DataType, typename FunctionType>
LookupTable<IndexType, DataType, FunctionType>::~LookupTable() 
{
    delete [] mp_indices;
    delete [] mp_data;
    
    mp_indices = NULL;
    mp_data    = NULL;
}

//==============================================================================

template<typename IndexType, typename DataType, typename FunctionType>
void LookupTable<IndexType, DataType, FunctionType>::initialize(
    int nrows, int nfuncs)
{
    m_num_indices   = nrows;
    m_num_functions = nfuncs;
    
    m_row_size = sizeof(DataType) * m_num_functions;
    
    m_is_constant_delta = false;
    
    // Allocate table storage
    mp_indices = new IndexType [m_num_indices];
    mp_data    = new DataType  [m_num_indices * m_num_functions];
    
    // Set initial values of all data to zero
    for (int i = 0; i < m_num_indices; ++i)
        mp_indices[i] = static_cast<IndexType>(0);
    
    for (int i = 0; i < m_num_indices * m_num_functions; ++i)
        mp_data[i] = static_cast<DataType>(0);
} // initialize()

//==============================================================================

template<typename IndexType, typename DataType, typename FunctionType>
void LookupTable<IndexType, DataType, FunctionType>::loadTable(
    const std::string &file_name) 
{
    // Open file for reading
    std::ifstream file(file_name.c_str(), std::ios::in);
    
    if (!file) throw FileNotFoundError(file_name);
    
    // Read the table dimensions initialize table data
    int num_indices;
    int num_functions;
    
    file >> num_indices >> num_functions;
    initialize(num_indices, num_functions);    
    
    m_is_constant_delta = true;
    
    // Read in each row of the table
    for (unsigned int i = 0; i < m_num_indices; i++) {
        // Read index of row i
        file >> mp_indices[i];
        
        // Ensure that indices are ordered in ascending order and do not repeat
        if (i > 0 && mp_indices[i] <= mp_indices[i-1]) {
            throw FileParseError(file_name, i+1)
                << "Lookup table is not ordered!\n"
                << "Index of row " << i+1 << " is <= to row " << i << ".";
        }
        
        // Test to see if the change in indices remains constant each row
        if (i > 1) {
            m_is_constant_delta &= 
                (abs((mp_indices[i]   - mp_indices[i-1]) - 
                     (mp_indices[i-1] - mp_indices[i-2])) < 
                 1.0e-6 * mp_indices[i-2]);
        }
        
        // Read in the function values associated with index i
        for (unsigned int j = 0; j < m_num_functions; j++)
            file >> mp_data[i*m_num_functions + j];
    }
    
    // Close the file
    file.close();
    
    // Write simple message showing that the table is loaded and give its
    // properties
    std::cout << std::endl;
    std::cout << "Loaded lookup table " + file_name + "." << std::endl;
    std::cout << "  min index      ~ " << minIndex() << std::endl;
    std::cout << "  max index      ~ " << maxIndex() << std::endl;
    std::cout << "  constant delta ~ " << 
        (m_is_constant_delta ? "True" : "False") << std::endl;
    std::cout << std::endl;
} // loadTable()

//==============================================================================

template<typename IndexType, typename DataType, typename FunctionType>
void LookupTable<IndexType, DataType, FunctionType>::populateTable(
    std::list< std::pair<IndexType, DataType *> > &table,
    typename std::list< std::pair<IndexType, DataType *> >::iterator &high,
    int nfuncs, double error_tolerance, InterpolationScheme scheme) const
{
    const int num_points = 3;
    
    typename std::list< std::pair<IndexType, DataType *> >::iterator low = high;
    low--;
    
    IndexType delta = (high->first - low->first) / (num_points + 1);
    
    // Loop over potential data points and determine index with maximum error
    // between interpolated and exact function values
    IndexType index;
    IndexType max_error_index;
    
    double ratio;
    double error;
    double max_error = 0.0;
    
    DataType *p_exact_values  = new DataType [nfuncs];
    DataType *p_interp_values = new DataType [nfuncs];
    
    FunctionType function;
    
    for (int i = 1; i <= num_points; ++i) {
        // Compute index of this point
        index = low->first + delta * static_cast<IndexType>(i);
        
        // Compute exact values
        //tableFunction(index, p_exact_values);
        function(index, p_exact_values, provider());
        
        // Now interpolated values
        ratio = static_cast<double>(index - low->first) / 
                static_cast<double>(high->first - low->first);
        
        switch(scheme) {
            case NEAREST_INDEX:
                _interpolateNearestIndex(
                    ratio, low->second, high->second, p_interp_values, 0, nfuncs,
                    Mutation::Numerics::Equals<DataType>());
                break;
            case LINEAR:
                _interpolateLinear(
                    ratio, low->second, high->second, p_interp_values, 0, nfuncs,
                    Mutation::Numerics::Equals<DataType>());
                break;
            case EXPONENTIAL:
                _interpolateExponential(
                    ratio, low->second, high->second, p_interp_values, 0, nfuncs,
                    Mutation::Numerics::Equals<DataType>());
                break;
        }
        
        // Compute error between interpolated and exact values
        for (int j = 0; j < nfuncs; ++j) {
            if (std::abs(p_exact_values[j]) < static_cast<DataType>(1.0e-10))
                error = std::abs(p_interp_values[j]);
            else 
                error = std::abs(static_cast<double>(
                    p_interp_values[j]/p_exact_values[j]) - 1.0);
            
            if (error > max_error) {
                max_error = error;
                max_error_index = index;
            }
        }
    }
    
    // We can go ahead and free p_interp_values so that memory does not build up
    // as recursion progresses (note that we cannot free p_exact_values because
    // it may be used to add data to table)
    delete [] p_interp_values;
    
    // If the max error found in interior points is greater than tolerance then
    // add that point to the table and repeat for regions on either side of new
    // point
    if (max_error > error_tolerance) {
        index = max_error_index;
        
        function(index, p_exact_values, provider());
        table.insert(
            high, std::pair<IndexType, DataType *>(index, p_exact_values));
        
        typename std::list<std::pair<IndexType, DataType *> >::iterator mid = 
            high;
        mid--;

        populateTable(
            table, mid, nfuncs, error_tolerance, scheme);
        populateTable(
            table, high, nfuncs, error_tolerance, scheme);
    } else {
        delete [] p_exact_values;
    }
} // populateTable()

//==============================================================================

template<typename IndexType, typename DataType, typename FunctionType>
void LookupTable<IndexType, DataType, FunctionType>::save(
    const std::string &file_name) const
{
    // Open file for writing
    std::ofstream file(file_name.c_str(), std::ios::out);
    
    if (!file) throw FileNotFoundError(file_name);
    
    // Set file format flags
    file.setf(std::ios::scientific | std::ios::right);
    file.precision(7);
    
    // Write out the table size dimensions
    file << std::setw(15) << m_num_indices;
    file << std::setw(15) << m_num_functions << std::endl;
    
    // Now write out all of the data
    for (int i = 0; i < m_num_indices; ++i) {
        file << std::setw(15) << mp_indices[i];
        
        for (int j = 0; j < m_num_functions; ++j)
            file << std::setw(15) << mp_data[i * m_num_functions + j];
        
        file << std::endl;
    }
    
    file << std::endl;
    
    // Close the file
    file.close();
} // save()

//==============================================================================

template<typename IndexType, typename DataType, typename FunctionType>
template<typename OP>
void LookupTable<IndexType, DataType, FunctionType>::lookup(
    const IndexType &index, const int start, const int end,
    DataType * const p_values, const OP& op, InterpolationScheme scheme) const
{
    unsigned int lower_row = 0;
    unsigned int upper_row = 1;
        
    // Search for the indices that bound the lookup index
    if (index > minIndex()) {
        if (m_is_constant_delta) {
            // With constant delta we can use a simple hash function lookup
            lower_row = static_cast<unsigned int>(
                (index - minIndex()) / (mp_indices[1] - mp_indices[0]));
            
            if (lower_row >= m_num_indices - 1)
                lower_row = m_num_indices - 2;
            
            upper_row = lower_row + 1;
        } else {
            // Otherwise do binary search for worst case O(log(n)) comparisons
            IndexType* p_upper_index = 
                std::lower_bound(mp_indices, mp_indices + m_num_indices - 1, index);
            
            upper_row = static_cast<unsigned int>(p_upper_index - mp_indices);
            lower_row = upper_row - 1;
        }
    }

//    std::cout << "Lookup table: lookup()" << std::endl;
//    std::cout << "# of indices = " << m_num_indices << std::endl;
//    std::cout << "x = " << index << std::endl;
//    std::cout << "index[" << lower_row << "]= " << mp_indices[lower_row] << std::endl;
//    std::cout << "index[" << upper_row << "]= " << mp_indices[upper_row] << std::endl;
    
    // Now compute the values of the functions at index using the specified
    // interpolation scheme
    double ratio = static_cast<double>(index - mp_indices[lower_row]) / 
                   (mp_indices[upper_row] - mp_indices[lower_row]);
    
    DataType *p_y1 = mp_data + lower_row * m_num_functions;
    DataType *p_y2 = mp_data + upper_row * m_num_functions;
    
    switch (scheme) {    
        case NEAREST_INDEX:
            _interpolateNearestIndex(
                ratio, p_y1, p_y2, p_values, start, end, op);
            break;   
                 
        case LINEAR:
            _interpolateLinear(
                ratio, p_y1, p_y2, p_values, start, end, op);
            break;   
                 
        case EXPONENTIAL:
            _interpolateExponential(
                ratio, p_y1, p_y2, p_values, start, end, op);
            break;
    }
} // lookup()

    } // namespace Utilities
} // namespace Mutation

#endif // LOOKUPTABLE_H

