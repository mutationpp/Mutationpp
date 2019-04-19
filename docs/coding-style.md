<a id="top"></a>
# Coding Style Guidelines

**Contents** <br>
- [1. General Guidelines](#1-general-guidelines)
- [2. C++](#2-c)
    - [2.1 Header files](#21-header-files)
    - [2.2 Source files](#22-source-files)
    - [2.3 Splitting over multiple lines](#23-splitting-over-multiple-lines)
    - [2.4 Math and logic](#24-math-and-logic)
- [3. Fortran](#3-fortran)

## 1. General Guidelines
Please follow these general guidelines when working in __ANY__ programming language.
- Maximum 80 characters per line
- Normal indentation is 4 spaces for each logical level
- Use spaces for indentation, not tabs
- Avoid trailing whitespace
- All files should be documented according to the [documentation guidelines](documenting-code.md#top)
- Use only standard language features

## 2. C++
All C++ code should be placed in either a `header file` or `source file` according to the following guidelines.
- Classes 
    - Name classes using [CamelCase](https://en.wikipedia.org/wiki/CamelCase), first letter capitalized
    - Place class definitions in a corresponding `header file` using the same name with the following structure

    ```c++
    // Class description using Doxygen
    class ClassName : public BaseClass
    {
    public:
        
        // Public constructors, destructors
        // Getters and setters (inlined when possible)
        // Other public methods
    
    private:
    
        // Private methods / variables
    
    }; // ClassName
    ```

    - Class members should be implemented in a separate `source file` with the same name
    - Small helper classes and/or child classes may be included in the main/base class header and source file when appropriate, each class should be documented separately according to the [documentation guidelines](documenting-code.md#top)
    - Use `const` whenever possible
    - Avoid explicitly inlining functions
- Control structures (`if`, `else`, `for`, `while`, etc.)
    - Do not use brackets for single line control structures
    - Use single letter variable names for loop counters whenever possible, ie: `i`, `j`, ...
- Variables
    - Use sensible names with [snake_case](https://en.wikipedia.org/wiki/Snake_case)
    - Avoid excessively long variable names 
    - Prefix class member variables with with `m_`
    - Prefix pointer variables with `p_`, or `mp_` for member pointers
- Functions
    - Name functions using [CamelCase](https://en.wikipedia.org/wiki/CamelCase), first letter lowercase 
    - Use the following formatting for function definitions
    
    ```c++
    // Short function like a getter or setter
    void shortFunction() { /* function body */ }
    
    // Single-line function
    void singleLineFunction() {
        // function body
    }
    
    // Multi-line function 
    void longFunction() 
    {
        // function body
    }


### 2.1 Header files
- Use only the `.h` extension, not `.H` or `.hpp`
- All header files must use `include guards`, use the name of the file in all caps with underscores to split words, followed by "_H"

```c++
#ifndef NAME_OF_FILE_H
#define NAME_OF_FILE_H
// body of file
#endif // NAME_OF_FILE_H
```

- Avoid function bodies in the header file except for templates and short inline functions
- Place class member template function definitions outside of the class declaration
- Do not use `using namespace` in header files, write out the full namespace paths

### 2.2 Source files
- Use only the `.cpp` extension, not `.cxx`
- Separate function definitions with 

```c++
//==============================================================================
```

### 2.3 Splitting over multiple lines
- Split long function headers on open parenthesis, indent following lines
- Do not put `const` on its own line

```c++
const Mutation::LongReturnType& Mutation::ClassName::longFunctionName(
    const VarType& long_var_name1, const VarType& long_var_name2) const 
{
    // function body
}
```

- Next split after return type if still too long, do not indent the function name

```c++
const Mutation::LongReturnType& 
Mutation::ClassName::longFunctionName(
    const VarType& long_var_name1, 
    const VarType& long_var_name2) const 
{
    // function body
}
```

- Split long formulae at `=` and operators, indent additional lines by 4 spaces

```c++
variable_name =
    a*(b + c) - c*(a + b)
```

- Split long `if` statements at logical operators, align left

```c++
if ((a == b) || (a == c) &&
    (b != c)) {
// if body
}
```

### 2.4 Math and logic
- Operator spacing

```c++
a + b, a - b
a*b, a/b
a & b, a ^ b
a = b, a != b
a < b, a > b, a >= b, a <= b
a || b, a && b
```

## 3. Fortran
- All program units should contain `implicit none` so that all variables, parameters, and functions must be declared
- Explicitly list module features when using a module

```fortran
use ModuleName, only : var1, var2
```
