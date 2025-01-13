# ifndef INTERPOLATION_H_
# define INTERPOLATION_H_

namespace interp
{
//----------------------------------------//
//------------Common Functions------------//
//--------Used by All Interpolators-------//
//----------------------------------------//
int find_segment(double, double*, int, int &);  // Find k such that x[k] <= x <= x[k+1]
double in_interval(double, double, double);     // Return x within x_min <= x <= x_max

//-------------------------------------//
//-----------One-Dimensional-----------//
//---------Linear Interpolation--------//
//-------------------------------------//
struct linear_spline_1D{
    int length;            // Length of input files (x and f(x))
    int accel;          // Index of previous table look up; used to increase spline speed
    double* x_vals;     // 1D array of x values
    double* f_vals;     // 1D array of f(x) values, f_vals[i] = f(x[i])
    double* slopes;     // Slopes used to generate natural cubic spline solution
};

void prep(struct linear_spline_1D &, int);   // Build arrays for an interpolation to get started
void set(struct linear_spline_1D &);         // Calculate slopes to set interpolation for evaluation
void clear(struct linear_spline_1D &);       // Clear arrays when finished

double eval_f(double, struct linear_spline_1D &);    // Evaluate f(x)
double eval_df(double, struct linear_spline_1D &);   // Evaluate df/dx

//-------------------------------------//
//-----------One-Dimensional-----------//
//---------Cubic Interpolation---------//
//-------------------------------------//
struct natural_cubic_spline_1D{
    int length;			// Length of input files (x and f(x))
    int accel;          // Index of previous table look up; used to increase spline speed
    double* x_vals;     // 1D array of x values
    double* f_vals;     // 1D array of f(x) values, f_vals[i] = f(x[i])
    double* slopes;     // Slopes used to generate natural cubic spline solution
};

void prep(struct natural_cubic_spline_1D &, int);   // Build arrays for an interpolation to get started
void set(struct natural_cubic_spline_1D &);         // Calculate slopes to set interpolation for evaluation
void clear(struct natural_cubic_spline_1D &);       // Clear arrays when finished

double eval_f(double, struct natural_cubic_spline_1D &);    // Evaluate f(x)
double eval_df(double, struct natural_cubic_spline_1D &);   // Evaluate df/dx
double eval_ddf(double, struct natural_cubic_spline_1D &);  // Evaluate d^2f/dx^2
double eval_dddf(double, struct natural_cubic_spline_1D &);  // Evaluate d^3f/dx^3

void eval_all(double, struct natural_cubic_spline_1D &, double &, double &);                                // Evaluate f(x) and f'(x)
void eval_all(double, struct natural_cubic_spline_1D &, double &, double &, double &);                      // Evaluate f(x), f'(x), and f''(x)
void eval_all(double, struct natural_cubic_spline_1D &, double &, double &, double &, double & );           // Evaluate f(x), f'(x), f''(x), and f'''(x)

//-------------------------------------//
//-----------Two-Dimensional-----------//
//------Multi-Stage Interpolation------//
//-------------------------------------//
struct natural_cubic_spline_2D{
    int length_x;           // Number of x nodes
    int length_y;           // Number of y nodes
    int accel[2];          // Indices of x,y point
    
    double* x_vals;         // 1D array of x values
    double* y_vals;         // 1D array of y values
    
    double** f_vals;       // 2D array of f(x,y) values, f_vals[i][j] = f(x[i], y[j])
    double** f_slopes;     // Slopes used to generate natural cubic spline solution at each x = constant node, describes f(x[i], y)
};

void prep(struct natural_cubic_spline_2D &, int, int);  // Build arrays for an interpolation to get started
void set(struct natural_cubic_spline_2D &);             // Calculate slopes to set interpolation for evaluation
void clear(struct natural_cubic_spline_2D &);           // Clear arrays when finished

double eval_node_f(double, struct natural_cubic_spline_2D, int, int);               // Evaluate the vertical spline of f at x_i
double eval_node_dfdy(double, struct natural_cubic_spline_2D, int, int);            // Evaluate the vertical spline of df/dy at x_i
double eval_node_ddfdydy(double, struct natural_cubic_spline_2D, int, int);         // Evaluate the vertical spline of d^2f/dy^2 at x_i
double eval_node_dddfdydydy(double, struct natural_cubic_spline_2D, int, int);      // Evaluate the vertical spline of d^3f/dy^3 at x_i

void set_node_slopes(double [], double [], struct natural_cubic_spline_2D);         // Set slopes along evaluated nodes

double eval_f(double, double, struct natural_cubic_spline_2D &);                    // Evaluate the spline at x, y
double eval_df(double, double, int, struct natural_cubic_spline_2D &);              // Evaluate df/dx_i at x, y
double eval_ddf(double, double, int, int, struct natural_cubic_spline_2D &);        // Evaluate d^2 f/dx_i dx_j at x, y 
double eval_dddf(double, double, int, int, int, struct natural_cubic_spline_2D &);  // Evaluate d^3 f/dx_i dx_j dx_k at x, y

void eval_all(double, double, struct natural_cubic_spline_2D &, double &, double []);
void eval_all(double, double, struct natural_cubic_spline_2D &, double &, double [], double []);
void eval_all(double, double, struct natural_cubic_spline_2D &, double &, double [], double [], double []);

//-------------------------------------//
//-----------Two-Dimensional-----------//
//--------Bicubic Interpolation--------//
//-------------------------------------//
extern double bicubic_conversion_matrix [16][16];

struct bicubic_spline{
    int length_x;       // Number of x nodes
    int length_y;       // Number of y nodes
    int accel[2];       // Indices of x,y point for acceleration
    
    double* x_vals;     // 1D array of x values
    double* y_vals;     // 1D array of y values
    
    double** f_vals;    // 2D array of f(x,y) values, f_vals[i][j] = f(x[i], y[j])
    double** dfdx_vals; // 2D array of df/dx(x,y) values, fdx_vals[i][j] = df/dy(x[i], y[j])
    double** dfdy_vals; // 2D array of df/dy(x,y) values, fdy_vals[i][j] = df/dy(x[i], y[j])
    
    double f_coeffs[4][4];      // Array of coefficients to compute f
    double dfdx_coeffs[4][4];   // Array of coefficients to compute df/dx
    double dfdy_coeffs[4][4];   // Array of coefficients to compute df/dy
};

void prep(struct bicubic_spline &, int, int);  // Build arrays for an interpolation to get started
void set(struct bicubic_spline &);             // Calculate slopes to set interpolation for evaluation
void clear(struct bicubic_spline &);           // Clear arrays when finished

bool same_cell(double, double, struct bicubic_spline);                    // Check if x,y are in the current acceleration cell ids
void update_spline_coeffs(int, int, struct bicubic_spline & );            // Update coefficients when moving between grid cells

double eval_f(double, double, struct bicubic_spline &);             // Evaluate the spline at x, y
double eval_df(double, double, int, struct bicubic_spline &);       // Evaluate df/dx_i at x, y
double eval_ddf(double, double, int, int, struct bicubic_spline &); // Evaluate d^2 f/dx_i dx_j at x, y

//---------------------------------------//
//-----------Three-Dimensional-----------//
//----------Hybrid Interpolation---------//
//---------------------------------------//
struct hybrid_spline_3D{
	int length_x;           // Number of x nodes
    int length_y;           // Number of y nodes
    int length_z;           // Number of z nodes
    int accel[3];           // Indices of x,y,z point

	double* x_vals;         // 1D array of x values
    double* y_vals;         // 1D array of y values
    double* z_vals;         // 1D array of z values
	
    double*** f_vals;       // 3D array of f(x,y,z) values, f_vals[i][j][k] = f(x[i], y[j], z[k])
	double*** f_slopes;     // Slopes used to generate natural cubic spline solution at each x and y = constant node, describes f(x[i], y[i], z)
    double*** dfdx_slopes;  // Slopes used to generate df/dx at each node point, describes df/dx @ x[i], y[j], z
    double*** dfdy_slopes;  // Slopes used to generate df/dy at each node point, describes df/dx @ x[i], y[j], z
};

void prep(struct hybrid_spline_3D &, int, int, int); // Build arrays for an interpolation to get started
void set(struct hybrid_spline_3D &);                 // Calculate slopes to set interpolation for evaluation
void clear(struct hybrid_spline_3D &);               // Clear arrays when finished

double eval_node_f(double, struct hybrid_spline_3D, int, int, int);         // Evaluate the vertical spline of f at x_i, y_j
double eval_node_dfdx(double, struct hybrid_spline_3D, int, int, int);      // Evaluate the vertical spline of df/dx at x_i, y_j
double eval_node_dfdy(double, struct hybrid_spline_3D, int, int, int);      // Evaluate the vertical spline of df/dy at x_i, y_j
double eval_node_dfdz(double, struct hybrid_spline_3D, int, int, int);      // Evaluate the vertical spline of df/dz at x_i, y_j
double eval_node_ddfdxdz(double, struct hybrid_spline_3D, int, int, int);   // Evaluate the vertical spline of d^2f/dxdz at x_i, y_j
double eval_node_ddfdydz(double, struct hybrid_spline_3D, int, int, int);   // Evaluate the vertical spline of d^2f/dydz at x_i, y_j
double eval_node_ddfdzdz(double, struct hybrid_spline_3D, int, int, int);   // Evaluate the vertical spline of d^2f/dz^2 at x_i, y_j

double finite_diff_dfdx(double, struct hybrid_spline_3D, int, int, int);        // Evaluate df/dx using finite differences
double finite_diff_dfdy(double, struct hybrid_spline_3D, int, int, int);        // Evaluate df/dy using finite differences

double finite_diff_ddfdxdx(double, struct hybrid_spline_3D, int, int, int);     // Evaluate d^2f/dx^2 using finite differences
double finite_diff_ddfdydy(double, struct hybrid_spline_3D, int, int, int);     // Evaluate d^2f/dy^2 using finite differences
double finite_diff_ddfdxdy(double, struct hybrid_spline_3D, int, int, int);     // Evaluate d^2f/dxdy using finite differences
double finite_diff_ddfdxdz(double, struct hybrid_spline_3D, int, int, int);     // Evaluate d^2f/dxdz using finite differences
double finite_diff_ddfdydz(double, struct hybrid_spline_3D, int, int, int);     // Evaluate d^2f/dydz using finite differences
    
double finite_diff_dddfdxdxdy(double, struct hybrid_spline_3D, int, int, int);   // Evaluate d^3f/dx^2dy using discrete derivatives
double finite_diff_dddfdxdydy(double, struct hybrid_spline_3D, int, int, int);   // Evaluate d^3f/dxdy^2 using discrete derivatives
double finite_diff_dddfdxdydz(double, struct hybrid_spline_3D, int, int, int);   // Evaluate d^3f/dxdydz using discrete derivatives
double finite_diff_dddfdxdzdz(double, struct hybrid_spline_3D, int, int, int);   // Evaluate d^3f/dxdz^2 using discrete derivatives
double finite_diff_dddfdydzdz(double, struct hybrid_spline_3D, int, int, int);   // Evaluate d^3f/dydz^2 using discrete derivatives

double finite_diff_ddddfdxdydzdz(double, struct hybrid_spline_3D, int, int, int);    // Evaluate d^4f/dxdydz^2 using discrete derivatives

double eval_f(double, double, double, struct hybrid_spline_3D &);               // Evaluate the spline at x, y, z
double eval_df(double, double, double, int, struct hybrid_spline_3D &);         // Evaluate df/dx_i at x, y, z
double eval_ddf(double, double, double, int, int, struct hybrid_spline_3D &);   // Evaluate d^2 f/dx_i dx_j at x, y, z

// Evaluate f, df/dx, df/dy, df/dz, d^2f/dx^2, d^2fdy^2, d^2f/dz^2, d^2f/dxdy, d^2f/dxdz, and d^2f/dydz
void eval_all(double, double, double, struct hybrid_spline_3D &, double &, double []);
void eval_all(double, double, double, struct hybrid_spline_3D &, double &, double [], double []);

}
#endif /* INTERPOLATION_H_ */
