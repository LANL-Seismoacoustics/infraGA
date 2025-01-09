# ifndef _INTERPOLATION_CPP_
# define _INTERPOLATION_CPP_

#include <iostream>
#include <math.h>
#include "interpolation.h"

using namespace std;

//----------------------------------------//
//------------Common Functions------------//
//--------Used by All Interpolators-------//
//----------------------------------------//
int interp::find_segment(double x, double* x_vals, int length, int & prev){
    bool done = false;
    
    if(x > x_vals[length - 1] || x < x_vals[0]){
        cout << "Cannot interpolate outside of given bounds.  x = " << x;
        cout << " is outside of bounds (" << x_vals[0] << ", " << x_vals[length - 1] << ")." << '\n';
    } else {
        if(x >= x_vals[prev] && x <= x_vals[prev+1]){
            done = true;
        }
        if(!done && prev + 2 <= length - 1){
            if(x >= x_vals[prev + 1] && x <= x_vals[prev + 2]){
                prev++;
                done = true;
            }
        }
        if(!done && prev - 1 >= 0){
            if(x >= x_vals[prev - 1] && x <= x_vals[prev]){
                prev--;
                done = true;
            }
        }
        if(!done){
            for(int n = 0; n < length; n++){
                if(x >= x_vals[n] && x <= x_vals[n + 1]){
                    prev = n;
                    break;
                } else if(x >= x_vals[(length - 1) - (n + 1)] && x < x_vals[(length - 1) - n]){
                    prev = (length - 1) - (n + 1);
                    break;
                }
            }
        }
    }
    return prev;
}


double interp::in_interval(double x, double x_min, double x_max){
    return min(max(x, x_min), x_max);
}


//-------------------------------------//
//-----------One-Dimensional-----------//
//---------Linear Interpolation--------//
//-------------------------------------//
void interp::prep(struct linear_spline_1D & spline, int length){
    spline.length = length;    spline.accel = 0;
    spline.x_vals = new double [spline.length];
    spline.f_vals = new double [spline.length];
    spline.slopes = new double [spline.length];
}

void interp::set(struct linear_spline_1D & spline){
    for(int i = 0; i < spline.length - 1; i++) {
        spline.slopes[i] = (spline.f_vals[i + 1] - spline.f_vals[i]) / (spline.x_vals[i + 1] - spline.x_vals[i]);
    }
    spline.slopes[spline.length - 1] = spline.slopes[spline.length - 2];
}

void interp::clear(struct linear_spline_1D & spline){
    delete spline.x_vals;
    delete spline.f_vals;
    delete spline.slopes;
}

double interp::eval_f(double x, struct linear_spline_1D & spline){
    int k = find_segment(x, spline.x_vals, spline.length, spline.accel);
    return spline.f_vals[k] + (x - spline.x_vals[k]) * spline.slopes[k];
}

double interp::eval_df(double x, struct linear_spline_1D & spline){
    int k = find_segment(x, spline.x_vals, spline.length, spline.accel);
    return spline.slopes[k];
}


//-------------------------------------//
//-----------One-Dimensional-----------//
//---------Cubic Interpolation---------//
//-------------------------------------//
void interp::prep(struct natural_cubic_spline_1D & spline, int length){
    spline.length = length;    spline.accel = 0;
    spline.x_vals = new double [spline.length];
    spline.f_vals = new double [spline.length];
    spline.slopes = new double [spline.length];
}

void interp::set(struct natural_cubic_spline_1D & spline){
    double ai, bi, ci, di;
    double new_c[spline.length-1];
    double new_d[spline.length];
    
    bi = 2.0 / (spline.x_vals[1] - spline.x_vals[0]);
    ci = 1.0 / (spline.x_vals[1] - spline.x_vals[0]);
    di = 3.0 * (spline.f_vals[1] - spline.f_vals[0]) / pow(spline.x_vals[1] - spline.x_vals[0], 2);
    
    new_c[0] = ci / bi;
    new_d[0] = di / bi;
    
    for(int i = 1; i < spline.length - 1; i++) {
        double  dx1 = spline.x_vals[i] - spline.x_vals[i - 1],
                dx2 = spline.x_vals[i + 1] - spline.x_vals[i];
        
        ai = 1.0 / dx1; bi = 2.0 * (1.0 / dx1 + 1.0 / dx2); ci = 1.0 / dx2;
        di = 3.0 * ((spline.f_vals[i] - spline.f_vals[i - 1]) / pow(dx1, 2)
                  + (spline.f_vals[i + 1] - spline.f_vals[i]) / pow(dx2, 2));
        
        new_c[i] = ci / (bi - new_c[i - 1] * ai);
        new_d[i] = (di - new_d[i - 1] * ai) / (bi - new_c[i - 1] * ai);
    }
    
    ai = 1.0 / (spline.x_vals[spline.length - 1] - spline.x_vals[spline.length - 2]);
    bi = 2.0 / (spline.x_vals[spline.length - 1] - spline.x_vals[spline.length - 2]);
    di = 3.0 * (spline.f_vals[spline.length - 1] - spline.f_vals[spline.length - 2])
          / pow(spline.x_vals[spline.length - 1] - spline.x_vals[spline.length - 2], 2);
    
    new_d[spline.length - 1] = (di - new_d[spline.length - 2] * ai) / (bi - new_c[spline.length - 2] * ai);
    
    spline.slopes[spline.length - 1] = new_d[spline.length - 1];
    for(int i = spline.length - 2; i > -1; i--){ spline.slopes[i] = new_d[i] - new_c[i] * spline.slopes[i + 1];}
}

void interp::clear(struct natural_cubic_spline_1D & spline){
    delete spline.x_vals;
    delete spline.f_vals;
    delete spline.slopes;
}

double interp::eval_f(double x, struct natural_cubic_spline_1D & spline){
    int k;
    double dx, X, df, A, B;

    k = find_segment(x, spline.x_vals, spline.length, spline.accel);

    dx = spline.x_vals[k + 1] - spline.x_vals[k];
    X = (x - spline.x_vals[k]) / dx;

    df = spline.f_vals[k + 1] - spline.f_vals[k];
    A = spline.slopes[k] * dx - df;
    B = -spline.slopes[k+1] * dx + df;
        
    return spline.f_vals[k] + X * (df + (1.0 - X) * (A * (1.0 - X) + B * X));
}

double interp::eval_df(double x, struct natural_cubic_spline_1D & spline){
    int k;
    double dx, X, df, A, B;

    k = find_segment(x, spline.x_vals, spline.length, spline.accel);
    
    dx = spline.x_vals[k + 1] - spline.x_vals[k];
    X = (x - spline.x_vals[k]) / dx;

    df = spline.f_vals[k + 1] - spline.f_vals[k];
    A = spline.slopes[k] * dx - df;
    B = -spline.slopes[k+1] * dx + df;

    return df / dx + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X) / dx + X * (1.0 - X) * (B - A) / dx;
}

double interp::eval_ddf(double x, struct natural_cubic_spline_1D & spline){
    int k;
    double dx, X, df, A, B;

    k = find_segment(x, spline.x_vals, spline.length, spline.accel);
    
    dx = spline.x_vals[k + 1] - spline.x_vals[k];
    X = (x - spline.x_vals[k]) / dx;

    df = spline.f_vals[k + 1] - spline.f_vals[k];
    A = spline.slopes[k] * dx - df;
    B = -spline.slopes[k+1] * dx + df;

    return 2.0 * (B - 2.0 * A + (A - B) * 3.0 * X) / pow(dx, 2);
}

double interp::eval_dddf(double x, struct natural_cubic_spline_1D & spline){
    int k;
    double dx, X, df, A, B;

    k = find_segment(x, spline.x_vals, spline.length, spline.accel);
    
    dx = spline.x_vals[k + 1] - spline.x_vals[k];
    X = (x - spline.x_vals[k]) / dx;

    df = spline.f_vals[k + 1] - spline.f_vals[k];
    A = spline.slopes[k] * dx - df;
    B = -spline.slopes[k+1] * dx + df;

    return 6.0 * (A - B) / pow(dx, 3);
}

void interp::eval_all(double x, struct natural_cubic_spline_1D & spline, double & f, double & dfdx){
    int k;
    double dx, X, df, A, B;

    k = find_segment(x, spline.x_vals, spline.length, spline.accel);
    
    dx = spline.x_vals[k + 1] - spline.x_vals[k];
    X = (x - spline.x_vals[k]) / dx;

    df = spline.f_vals[k + 1] - spline.f_vals[k];
    A = spline.slopes[k] * dx - df;
    B = -spline.slopes[k+1] * dx + df;

    f = spline.f_vals[k] + X * (df + (1.0 - X) * (A * (1.0 - X) + B * X));
    dfdx = df / dx + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X) / dx + X * (1.0 - X) * (B - A) / dx;
}

void interp::eval_all(double x, struct natural_cubic_spline_1D & spline, double & f, double & dfdx, double & ddfdxdx){
    int k;
    double dx, X, df, A, B;

    k = find_segment(x, spline.x_vals, spline.length, spline.accel);
    
    dx = spline.x_vals[k + 1] - spline.x_vals[k];
    X = (x - spline.x_vals[k]) / dx;

    df = spline.f_vals[k + 1] - spline.f_vals[k];
    A = spline.slopes[k] * dx - df;
    B = -spline.slopes[k+1] * dx + df;

    f = spline.f_vals[k] + X * (df + (1.0 - X) * (A * (1.0 - X) + B * X));
    dfdx = df / dx + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X) / dx + X * (1.0 - X) * (B - A) / dx;
    ddfdxdx = 2.0 * (B - 2.0 * A + (A - B) * 3.0 * X) / pow(dx, 2);
}

void interp::eval_all(double x, struct natural_cubic_spline_1D & spline, double & f, double & dfdx, double & ddfdxdx, double & dddfdxdxdx){
    int k;
    double dx, X, df, A, B;

    k = find_segment(x, spline.x_vals, spline.length, spline.accel);
    
    dx = spline.x_vals[k + 1] - spline.x_vals[k];
    X = (x - spline.x_vals[k]) / dx;

    df = spline.f_vals[k + 1] - spline.f_vals[k];
    A = spline.slopes[k] * dx - df;
    B = -spline.slopes[k+1] * dx + df;

    f = spline.f_vals[k] + X * (df + (1.0 - X) * (A * (1.0 - X) + B * X));
    dfdx = df / dx + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X) / dx + X * (1.0 - X) * (B - A) / dx;
    ddfdxdx = 2.0 * (B - 2.0 * A + (A - B) * 3.0 * X) / pow(dx, 2);
    dddfdxdxdx = 6.0 * (A - B) / pow(dx, 3);
}

//-------------------------------------//
//-----------Two-Dimensional-----------//
//------Multi-Stage Interpolation------//
//-------------------------------------//
void interp::prep(struct natural_cubic_spline_2D & spline, int nx, int ny){
    spline.length_x = nx;   spline.accel[0] = 0;   spline.x_vals = new double [spline.length_x];
    spline.length_y = ny;   spline.accel[1] = 0;   spline.y_vals = new double [spline.length_y];
    
    spline.f_vals = new double* [spline.length_x];
    spline.f_slopes = new double* [spline.length_x];
    
    for(int nx = 0; nx < spline.length_x; nx++){
        spline.f_vals[nx] = new double [spline.length_y];
        spline.f_slopes[nx] = new double [spline.length_y];
    }
}

void interp::set(struct natural_cubic_spline_2D & spline){
    double dy1, dy2, ai, bi, ci, di;
    double new_c[spline.length_y - 1];
    double new_d[spline.length_y];
    
    for(int nx = 0; nx < spline.length_x; nx++){
        bi = 2.0 / (spline.y_vals[1] - spline.y_vals[0]);
        ci = 1.0 / (spline.y_vals[1] - spline.y_vals[0]);
        di = 3.0 * (spline.f_vals[nx][1] - spline.f_vals[nx][0]) / pow(spline.y_vals[1] - spline.y_vals[0], 2);
        
        new_c[0] = ci / bi;
        new_d[0] = di / bi;
        
        for(int ny = 1; ny < spline.length_y - 1; ny++) {
            dy1 = spline.y_vals[ny] - spline.y_vals[ny - 1],
            dy2 = spline.y_vals[ny + 1] - spline.y_vals[ny];
        
            ai = 1.0 / dy1; bi = 2.0 * (1.0 / dy1 + 1.0 / dy2); ci = 1.0 / dy2;
            di = 3.0 * ((spline.f_vals[nx][ny] - spline.f_vals[nx][ny - 1]) / pow(dy1, 2)
                      + (spline.f_vals[nx][ny + 1] - spline.f_vals[nx][ny]) / pow(dy2, 2));
            
            new_c[ny] = ci / (bi - new_c[ny - 1] * ai);
            new_d[ny] = (di - new_d[ny - 1] * ai) / (bi - new_c[ny - 1] * ai);
        }
        
        ai = 1.0 / (spline.y_vals[spline.length_y - 1] - spline.y_vals[spline.length_y - 2]);
        bi = 2.0 / (spline.y_vals[spline.length_y - 1] - spline.y_vals[spline.length_y - 2]);
        di = 3.0 * (spline.f_vals[nx][spline.length_y - 1] - spline.f_vals[nx][spline.length_y - 2])
                    / pow(spline.y_vals[spline.length_y - 1] - spline.y_vals[spline.length_y - 2], 2);
        
        new_d[spline.length_y - 1] = (di - new_d[spline.length_y - 2] * ai) / (bi - new_c[spline.length_y - 2] * ai);
        
        spline.f_slopes[nx][spline.length_y - 1] = new_d[spline.length_y - 1];
        for(int ny = spline.length_y - 2; ny >= 0; ny--){
            spline.f_slopes[nx][ny] = new_d[ny] - new_c[ny] * spline.f_slopes[nx][ny + 1];
        }


    }
}

void interp::clear(struct natural_cubic_spline_2D & spline){
    delete spline.x_vals;
    delete spline.y_vals;
    
    for(int nx = 0; nx < spline.length_x; nx++){
        delete spline.f_vals[nx];
        delete spline.f_slopes[nx];
    }
    delete spline.f_vals;
    delete spline.f_slopes;
}

double interp::eval_node_f(double y, struct natural_cubic_spline_2D spline, int nx, int ny){
    double df, dy, X, A, B;

    df = spline.f_vals[nx][ny + 1] - spline.f_vals[nx][ny];
    dy = spline.y_vals[ny + 1] - spline.y_vals[ny];
    
    X = (y - spline.y_vals[ny]) / dy;
    A = spline.f_slopes[nx][ny] * dy - df;
    B = -spline.f_slopes[nx][ny + 1] * dy + df;
    
    return spline.f_vals[nx][ny] + X * (df + (1.0 - X) * (A * (1.0 - X) + B * X));
}

double interp::eval_node_dfdy(double y, struct natural_cubic_spline_2D spline, int nx, int ny){
    double df, dy, X, A, B;

    df = spline.f_vals[nx][ny + 1] - spline.f_vals[nx][ny];
    dy = spline.y_vals[ny + 1] - spline.y_vals[ny];
    
    X = (y - spline.y_vals[ny]) / dy;
    A = spline.f_slopes[nx][ny] * dy - df;
    B = -spline.f_slopes[nx][ny + 1] * dy + df;

    return df / dy + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X) / dy + X * (1.0 - X) * (B - A) / dy;
}

double interp::eval_node_ddfdydy(double y, struct natural_cubic_spline_2D spline, int nx, int ny){
    double df, dy, X, A, B;

    df = spline.f_vals[nx][ny + 1] - spline.f_vals[nx][ny];
    dy = spline.y_vals[ny + 1] - spline.y_vals[ny];
    
    X = (y - spline.y_vals[ny]) / dy;
    A = spline.f_slopes[nx][ny] * dy - df;
    B = -spline.f_slopes[nx][ny + 1] * dy + df;

    return 2.0 * (B - 2.0 * A + (A - B) * 3.0 * X) / pow(dy, 2);
}

double interp::eval_node_dddfdydydy(double y, struct natural_cubic_spline_2D spline, int nx, int ny){
    double df, dy, X, A, B;

    df = spline.f_vals[nx][ny + 1] - spline.f_vals[nx][ny];
    dy = spline.y_vals[ny + 1] - spline.y_vals[ny];
    
    X = (y - spline.y_vals[ny]) / dy;
    A = spline.f_slopes[nx][ny] * dy - df;
    B = -spline.f_slopes[nx][ny + 1] * dy + df;

    return 6.0 * (A - B) / pow(dy, 3);
}

void interp::set_node_slopes(double node_vals [], double node_slopes [], struct natural_cubic_spline_2D spline){
    double dx1, dx2, ai, bi, ci, di;
    double new_c[spline.length_x - 1];
    double new_d[spline.length_x];

    bi = 2.0 / (spline.x_vals[1] - spline.x_vals[0]);
    ci = 1.0 / (spline.x_vals[1] - spline.x_vals[0]);
    di = 3.0 * (node_vals[1] - node_vals[0]) / pow(spline.x_vals[1] - spline.x_vals[0], 2);
    
    new_c[0] = ci / bi;
    new_d[0] = di / bi;
    
    for(int nx = 1; nx < spline.length_x - 1; nx++) {
        dx1 = spline.x_vals[nx] - spline.x_vals[nx - 1];
        dx2 = spline.x_vals[nx + 1] - spline.x_vals[nx];
        
        ai = 1.0 / dx1; bi = 2.0 * (1.0 / dx1 + 1.0 / dx2); ci = 1.0 / dx2;
        di = 3.0 * ((node_vals[nx] - node_vals[nx - 1]) / pow(dx1, 2)
                  + (node_vals[nx + 1] - node_vals[nx]) / pow(dx2, 2));
        
        new_c[nx] = ci / (bi - new_c[nx - 1] * ai);
        new_d[nx] = (di - new_d[nx - 1] * ai)/(bi - new_c[nx - 1] * ai);
    }
    
    ai = 1.0 / (spline.x_vals[spline.length_x - 1] - spline.x_vals[spline.length_x - 2]);
    bi = 2.0 / (spline.x_vals[spline.length_x - 1] - spline.x_vals[spline.length_x - 2]);
    di = 3.0 * (node_vals[spline.length_x - 1] - node_vals[spline.length_x - 2]) / pow(spline.x_vals[spline.length_x - 1] - spline.x_vals[spline.length_x - 2], 2);
    
    new_d[spline.length_x - 1] = (di - new_d[spline.length_x - 2] * ai) / (bi - new_c[spline.length_x - 2] * ai);

    node_slopes[spline.length_x - 1] = new_d[spline.length_x - 1];
    for(int nx = spline.length_x - 2; nx >= 0; nx--){
        node_slopes[nx] = new_d[nx] - new_c[nx] * node_slopes[nx + 1];
    }
}


double interp::eval_f(double x, double y, struct natural_cubic_spline_2D & spline){
    int kx, ky;
    double node_vals[spline.length_x], node_slopes[spline.length_x];
    double df, dx, X, A, B;

    kx = find_segment(x, spline.x_vals, spline.length_x, spline.accel[0]);
    ky = find_segment(y, spline.y_vals, spline.length_y, spline.accel[1]);

    for(int nx = 0; nx < spline.length_x; nx++){
        node_vals[nx] = eval_node_f(y, spline, nx, ky);
    }
    set_node_slopes(node_vals, node_slopes, spline);
    
    df = node_vals[kx + 1] - node_vals[kx];
    dx = spline.x_vals[kx + 1] - spline.x_vals[kx];
    
    X = (x - spline.x_vals[kx]) / dx;
    A = node_slopes[kx] * dx - df;
    B = -node_slopes[kx + 1] * dx + df;
    
    return node_vals[kx] + X * (df + (1.0 - X) * (A * (1.0 - X) + B * X));
}

double interp::eval_df(double x, double y, int n, struct natural_cubic_spline_2D & spline){
    int kx, ky;
    double node_vals[spline.length_x], node_slopes[spline.length_x];
    double df, dx, X, A, B;

    kx = find_segment(x, spline.x_vals, spline.length_x, spline.accel[0]);
    ky = find_segment(y, spline.y_vals, spline.length_y, spline.accel[1]);

    for(int nx = 0; nx < spline.length_x; nx++){ 
        if(n == 0){
            node_vals[nx] = eval_node_f(y, spline, nx, ky);
        } else {
            node_vals[nx] = eval_node_dfdy(y, spline, nx, ky);
        }
    }
    set_node_slopes(node_vals, node_slopes, spline);

    df = node_vals[kx + 1] - node_vals[kx];
    dx = spline.x_vals[kx + 1] - spline.x_vals[kx];

    X = (x - spline.x_vals[kx]) / dx;
    A = node_slopes[kx] * dx - df;
    B = -node_slopes[kx + 1] * dx + df;
    
    if(n == 0){
        return df / dx + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X) / dx + X * (1.0 - X) * (B - A) / dx;
    } else {
        return node_vals[kx] + X * (df + (1.0 - X) * (A * (1.0 - X) + B * X));
    }
}

double interp::eval_ddf(double x, double y, int n1, int n2, struct natural_cubic_spline_2D & spline){
    int kx, ky;
    double node_vals[spline.length_x], node_slopes[spline.length_x];
    double df, dx, X, A, B;

    kx = find_segment(x, spline.x_vals, spline.length_x, spline.accel[0]);
    ky = find_segment(y, spline.y_vals, spline.length_y, spline.accel[1]);
    
    for(int nx = 0; nx < spline.length_x; nx++){
        if(n1 + n2 == 0){
            node_vals[nx] = eval_node_f(y, spline, nx, ky);
        } else if(n1 + n2 == 1){
            node_vals[nx] = eval_node_dfdy(y, spline, nx, ky);
        } else {
            node_vals[nx] = eval_node_ddfdydy(y, spline, nx, ky);
        }
    }
    set_node_slopes(node_vals, node_slopes, spline);

    df = node_vals[kx + 1] - node_vals[kx];
    dx = spline.x_vals[kx + 1] - spline.x_vals[kx];

    X = (x - spline.x_vals[kx]) / dx;
    A = node_slopes[kx] * dx - df;
    B = -node_slopes[kx + 1] * dx + df;

    if(n1 + n2 == 0){
        return 2.0 * (B - 2.0 * A + (A - B) * 3.0 * X) / pow(dx, 2);
    } else if(n1  + n2 == 1){
        return df / dx + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X) / dx + X * (1.0 - X) * (B - A) / dx;
    } else { 
        return node_vals[kx] + X * (df + (1.0 - X) * (A * (1.0 - X ) + B * X));
    }
}

double interp::eval_dddf(double x, double y, int n1, int n2, int n3, struct natural_cubic_spline_2D & spline){
    int kx, ky;
    double node_vals[spline.length_x], node_slopes[spline.length_x];
    double df, dx, X, A, B;

    kx = find_segment(x, spline.x_vals, spline.length_x, spline.accel[0]);
    ky = find_segment(y, spline.y_vals, spline.length_y, spline.accel[1]);    

    for(int nx = 0; nx < spline.length_x; nx++){ 
        if(n1 + n2 + n3 == 0){
            node_vals[nx] = eval_node_f(y, spline, nx, ky);
        } else if(n1 + n2 + n3 == 1){
            node_vals[nx] = eval_node_dfdy(y, spline, nx, ky);
        } else if(n1 + n2 + n3 == 2){
            node_vals[nx] = eval_node_ddfdydy(y, spline, nx, ky);
        } else {
            node_vals[nx] = eval_node_dddfdydydy(y, spline, nx, ky);
        }
    }
    set_node_slopes(node_vals, node_slopes, spline);
    df = node_vals[kx + 1] - node_vals[kx];
    dx = spline.x_vals[kx + 1] - spline.x_vals[kx];

    X = (x - spline.x_vals[kx]) / dx;
    A = node_slopes[kx] * dx - df;
    B = -node_slopes[kx + 1] * dx + df;

    if(n1 + n2 + n3 == 0){
        return 6.0 * (A - B) / pow(dx, 3);
    } else if(n1 + n2 + n3 == 1){
        return 2.0 * (B - 2.0 * A + (A - B) * 3.0 * X) / pow(dx, 2);
    } else if(n1 + n2 + n3 == 2){
        return df / dx + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X) / dx + X * (1.0 - X) * (B - A) / dx;
    } else {
        return node_vals[kx] + X * (df + (1.0 - X) * (A * (1.0 - X ) + B * X));
    }
}

void interp::eval_all(double x, double y, struct natural_cubic_spline_2D & spline, double & f, double dfdx []){
    int kx, ky;
    double node_f_vals[spline.length_x], node_f_slopes[spline.length_x];
    double node_df_vals[spline.length_x], node_df_slopes[spline.length_x];
    double df, dx, X, A, B;

    kx = find_segment(x, spline.x_vals, spline.length_x, spline.accel[0]);
    ky = find_segment(y, spline.y_vals, spline.length_y, spline.accel[1]);    

    for(int nx = 0; nx < spline.length_x; nx++){ 
        node_f_vals[nx] = eval_node_f(y, spline, nx, ky);
        node_df_vals[nx] = eval_node_dfdy(y, spline, nx, ky);
    }
    set_node_slopes(node_f_vals, node_f_slopes, spline);
    set_node_slopes(node_df_vals, node_df_slopes, spline);

    dx = spline.x_vals[kx + 1] - spline.x_vals[kx];
    X = (x - spline.x_vals[kx]) / dx;

    df = node_f_vals[kx + 1] - node_f_vals[kx];
    A = node_f_slopes[kx] * dx - df;
    B = -node_f_slopes[kx + 1] * dx + df;

    f = node_f_vals[kx] + X * (df + (1.0 - X) * (A * (1.0 - X) + B * X));
    dfdx[0] = df / dx + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X) / dx + X * (1.0 - X) * (B - A) / dx;

    df = node_df_vals[kx + 1] - node_df_vals[kx];
    A = node_df_slopes[kx] * dx - df;
    B = -node_df_slopes[kx + 1] * dx + df;

    dfdx[1] = node_df_vals[kx] + X * (df + (1.0 - X) * (A * (1.0 - X) + B * X));
}

void interp::eval_all(double x, double y, struct natural_cubic_spline_2D & spline, double & f, double dfdx [], double ddfdxdx []){
    int kx, ky;
    double node_f_vals[spline.length_x], node_f_slopes[spline.length_x];
    double node_df_vals[spline.length_x], node_df_slopes[spline.length_x];
    double node_ddf_vals[spline.length_x], node_ddf_slopes[spline.length_x];
    double df, dx, X, A, B;

    kx = find_segment(x, spline.x_vals, spline.length_x, spline.accel[0]);
    ky = find_segment(y, spline.y_vals, spline.length_y, spline.accel[1]);    

    for(int nx = 0; nx < spline.length_x; nx++){ 
        node_f_vals[nx] = eval_node_f(y, spline, nx, ky);
        node_df_vals[nx] = eval_node_dfdy(y, spline, nx, ky);
        node_ddf_vals[nx] = eval_node_ddfdydy(y, spline, nx, ky);
    }
    set_node_slopes(node_f_vals, node_f_slopes, spline);
    set_node_slopes(node_df_vals, node_df_slopes, spline);
    set_node_slopes(node_ddf_vals, node_ddf_slopes, spline);

    dx = spline.x_vals[kx + 1] - spline.x_vals[kx];
    X = (x - spline.x_vals[kx]) / dx;

    df = node_f_vals[kx + 1] - node_f_vals[kx];
    A = node_f_slopes[kx] * dx - df;
    B = -node_f_slopes[kx + 1] * dx + df;

    f = node_f_vals[kx] + X * (df + (1.0 - X) * (A * (1.0 - X) + B * X));
    dfdx[0] = df / dx + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X) / dx + X * (1.0 - X) * (B - A) / dx;
    ddfdxdx[0] = 2.0 * (B - 2.0 * A + (A - B) * 3.0 * X) / pow(dx, 2);

    df = node_df_vals[kx + 1] - node_df_vals[kx];
    A = node_df_slopes[kx] * dx - df;
    B = -node_df_slopes[kx + 1] * dx + df;

    dfdx[1] = node_df_vals[kx] + X * (df + (1.0 - X) * (A * (1.0 - X) + B * X));
    ddfdxdx[1] = df / dx + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X) / dx + X * (1.0 - X) * (B - A) / dx;

    df = node_ddf_vals[kx + 1] - node_ddf_vals[kx];
    A = node_ddf_slopes[kx] * dx - df;
    B = -node_ddf_slopes[kx + 1] * dx + df;

    ddfdxdx[2] = node_ddf_vals[kx] + X * (df + (1.0 - X) * (A * (1.0 - X) + B * X));
}


void interp::eval_all(double x, double y, struct natural_cubic_spline_2D & spline, double & f, double dfdx [], double ddfdxdx [], double dddfdxdxdx []){
    int kx, ky;
    double node_f_vals[spline.length_x], node_f_slopes[spline.length_x];
    double node_df_vals[spline.length_x], node_df_slopes[spline.length_x];
    double node_ddf_vals[spline.length_x], node_ddf_slopes[spline.length_x];
    double node_dddf_vals[spline.length_x], node_dddf_slopes[spline.length_x];
    double df, dx, X, A, B;

    kx = find_segment(x, spline.x_vals, spline.length_x, spline.accel[0]);
    ky = find_segment(y, spline.y_vals, spline.length_y, spline.accel[1]);    

    for(int nx = 0; nx < spline.length_x; nx++){ 
        node_f_vals[nx] = eval_node_f(y, spline, nx, ky);
        node_df_vals[nx] = eval_node_dfdy(y, spline, nx, ky);
        node_ddf_vals[nx] = eval_node_ddfdydy(y, spline, nx, ky);
        node_dddf_vals[nx] = eval_node_dddfdydydy(y, spline, nx, ky);
    }
    set_node_slopes(node_f_vals, node_f_slopes, spline);
    set_node_slopes(node_df_vals, node_df_slopes, spline);
    set_node_slopes(node_ddf_vals, node_ddf_slopes, spline);
    set_node_slopes(node_dddf_vals, node_dddf_slopes, spline);

    dx = spline.x_vals[kx + 1] - spline.x_vals[kx];
    X = (x - spline.x_vals[kx]) / dx;

    df = node_f_vals[kx + 1] - node_f_vals[kx];
    A = node_f_slopes[kx] * dx - df;
    B = -node_f_slopes[kx + 1] * dx + df;
    
    f = node_f_vals[kx] + X * (df + (1.0 - X) * (A * (1.0 - X) + B * X));
    dfdx[0] = df / dx + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X) / dx + X * (1.0 - X) * (B - A) / dx;
    ddfdxdx[0] = 2.0 * (B - 2.0 * A + (A - B) * 3.0 * X) / pow(dx, 2);
    dddfdxdxdx[0] = 6.0 * (A - B) / pow(dx, 3);

    df = node_df_vals[kx + 1] - node_df_vals[kx];
    A = node_df_slopes[kx] * dx - df;
    B = -node_df_slopes[kx + 1] * dx + df;

    dfdx[1] = node_df_vals[kx] + X * (df + (1.0 - X) * (A * (1.0 - X) + B * X));
    ddfdxdx[1] = df / dx + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X) / dx + X * (1.0 - X) * (B - A) / dx;
    dddfdxdxdx[1] = 2.0 * (B - 2.0 * A + (A - B) * 3.0 * X) / pow(dx, 2);

    df = node_ddf_vals[kx + 1] - node_ddf_vals[kx];
    A = node_ddf_slopes[kx] * dx - df;
    B = -node_ddf_slopes[kx + 1] * dx + df;

    ddfdxdx[2] = node_ddf_vals[kx] + X * (df + (1.0 - X) * (A * (1.0 - X) + B * X));
    dddfdxdxdx[2] = df / dx + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X) / dx + X * (1.0 - X) * (B - A) / dx;

    df = node_dddf_vals[kx + 1] - node_dddf_vals[kx];
    A = node_dddf_slopes[kx] * dx - df;
    B = -node_dddf_slopes[kx + 1] * dx + df;

    dddfdxdxdx[3] = node_dddf_vals[kx] + X * (df + (1.0 - X) * (A * (1.0 - X) + B * X));
}

//-------------------------------------//
//-----------Two-Dimensional-----------//
//--------Bicubic Interpolation--------//
//-------------------------------------//
double interp::bicubic_conversion_matrix [16][16] = {
    { 1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
    { 0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
    {-3.0,  3.0,  0.0,  0.0, -2.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
    { 2.0, -2.0,  0.0,  0.0,  1.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0},
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -3.0,  3.0,  0.0,  0.0, -2.0, -1.0,  0.0,  0.0},
    { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  2.0, -2.0,  0.0,  0.0,  1.0,  1.0,  0.0,  0.0},
    {-3.0,  0.0,  3.0,  0.0,  0.0,  0.0,  0.0,  0.0, -2.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0},
    { 0.0,  0.0,  0.0,  0.0, -3.0,  0.0,  3.0,  0.0,  0.0,  0.0,  0.0,  0.0, -2.0,  0.0, -1.0,  0.0},
    { 9.0, -9.0, -9.0,  9.0,  6.0,  3.0, -6.0, -3.0,  6.0, -6.0,  3.0, -3.0,  4.0,  2.0,  2.0,  1.0},
    {-6.0,  6.0,  6.0, -6.0, -3.0, -3.0,  3.0,  3.0, -4.0,  4.0, -2.0,  2.0, -2.0, -2.0, -1.0, -1.0},
    { 2.0,  0.0, -2.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0},
    { 0.0,  0.0,  0.0,  0.0,  2.0,  0.0, -2.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  1.0,  0.0},
    {-6.0,  6.0,  6.0, -6.0, -4.0, -2.0,  4.0,  2.0, -3.0,  3.0, -3.0,  3.0, -2.0, -1.0, -2.0, -1.0},
    { 4.0, -4.0, -4.0,  4.0,  2.0,  2.0, -2.0, -2.0,  2.0, -2.0,  2.0, -2.0,  1.0,  1.0,  1.0,  1.0}
};

void interp::prep(struct bicubic_spline & spline, int nx, int ny){
    spline.length_x = nx;
    spline.length_y = ny;

    spline.accel[0] = 0;
    spline.accel[1] = 0;

    spline.x_vals = new double [spline.length_x];
    spline.y_vals = new double [spline.length_y];

    spline.f_vals = new double* [spline.length_x];
    for(int i = 0; i < spline.length_x; i++){
        spline.f_vals[i] = new double [spline.length_y];
    }
}

void interp::set(struct bicubic_spline & spline){}

void interp::clear(struct bicubic_spline & spline){
    delete spline.x_vals;
    delete spline.y_vals;
    
    for(int i = 0; i < spline.length_x; i++){
        delete spline.f_vals[i];
    }
    delete spline.f_vals;
}

bool interp::same_cell(double x, double y, struct bicubic_spline spline){
    if(spline.x_vals[spline.accel[0]] <= x && x <= spline.x_vals[spline.accel[0] + 1]){
        if(spline.y_vals[spline.accel[1]] <= y && y <= spline.y_vals[spline.accel[1] + 1]){
            return true;
        }
    }
    return false;
}

void interp::update_spline_coeffs(int nx1, int ny1, struct bicubic_spline & spline){
    double A[16], X[16];
    int nx1_up, nx1_dn, nx2, nx2_up, nx2_dn, ny1_up, ny1_dn, ny2, ny2_up, ny2_dn;
    
    nx1_up = nx1 + 1; nx1_dn = max(nx1 - 1, 0); nx2 = nx1 + 1;
    ny1_up = ny1 + 1; ny1_dn = max(ny1 - 1, 0); ny2 = ny1 + 1;
    
    nx2_up = min(nx2 + 1, spline.length_x - 1); nx2_dn = nx2 - 1;
    ny2_up = min(ny2 + 1, spline.length_y - 1); ny2_dn = ny2 - 1;
    
    X[0] = spline.f_vals[nx1][ny1];                                         X[1] = spline.f_vals[nx2][ny1];
    X[2] = spline.f_vals[nx1][ny2];                                         X[3] = spline.f_vals[nx2][ny2];
    
    X[4] = spline.f_vals[nx1_up][ny1] - spline.f_vals[nx1_dn][ny1];         X[5] = spline.f_vals[nx2_up][ny1] - spline.f_vals[nx2_dn][ny1];
    X[6] = spline.f_vals[nx1_up][ny2] - spline.f_vals[nx1_dn][ny2];         X[7] = spline.f_vals[nx2_up][ny2] - spline.f_vals[nx2_dn][ny2];
    
    X[8] =  spline.f_vals[nx1][ny1_up] - spline.f_vals[nx1][ny1_dn];        X[9] =  spline.f_vals[nx2][ny1_up] - spline.f_vals[nx2][ny1_dn];
    X[10] = spline.f_vals[nx1][ny2_up] - spline.f_vals[nx1][ny2_dn];        X[11] = spline.f_vals[nx2][ny2_up] - spline.f_vals[nx2][ny2_dn];
    
    X[12] = spline.f_vals[nx1_up][ny1_up] - spline.f_vals[nx1_dn][ny1_dn];  X[13] = spline.f_vals[nx2_up][ny1_up] - spline.f_vals[nx2_dn][ny1_dn];
    X[14] = spline.f_vals[nx1_up][ny2_up] - spline.f_vals[nx1_dn][ny2_dn];  X[15] = spline.f_vals[nx2_up][ny2_up] - spline.f_vals[nx2_dn][ny2_dn];
    
    for(int j = 0; j < 16; j++){
        A[j] = 0.0;
        for(int k = 0; k < 16; k++){
            A[j] += bicubic_conversion_matrix[j][k] * X[k];
        }
    }
    for(int j= 0 ; j < 4; j++){
        for(int k = 0; k < 4; k++){
            spline.f_coeffs[j][k] = A[k * 4 + j];
        }
    }
}

double interp::eval_f(double x, double y, struct bicubic_spline & spline){
    int nx, ny;
    double x_scaled, y_scaled, result;
    
    if(!same_cell(x, y, spline)){
        nx = find_segment(x, spline.x_vals, spline.length_x, spline.accel[0]);
        ny = find_segment(y, spline.y_vals, spline.length_y, spline.accel[1]);
        update_spline_coeffs(nx, ny, spline);
    } else {
        nx = spline.accel[0];
        ny = spline.accel[1];
    }
    
    x_scaled = (x - spline.x_vals[nx]) / (spline.x_vals[nx + 1] - spline.x_vals[nx]);
    y_scaled = (y - spline.y_vals[ny]) / (spline.y_vals[ny + 1] - spline.y_vals[ny]);

    result = 0.0;
    for(int j = 0; j < 4; j++){
        for(int k = 0; k < 4; k++){
            result += pow(x_scaled, j) * pow(y_scaled, k) * spline.f_coeffs[j][k];
        }
    }
    return result;
}

double interp::eval_df(double x, double y, int n, struct bicubic_spline & spline){
    int nx, ny;
    double dx, dy, x_scaled, y_scaled, result = 0.0;
    
    if(!same_cell(x, y, spline)){
        nx = find_segment(x, spline.x_vals, spline.length_x, spline.accel[0]);
        ny = find_segment(y, spline.y_vals, spline.length_y, spline.accel[1]);
        update_spline_coeffs(nx, ny, spline);
    } else {
        nx = spline.accel[0];
        ny = spline.accel[1];
    }
    
    dx = spline.x_vals[nx + 1] - spline.x_vals[nx];     x_scaled = (x - spline.x_vals[nx]) / dx;
    dy = spline.y_vals[ny + 1] - spline.y_vals[ny];     y_scaled = (y - spline.y_vals[ny]) / dy;

    result = 0.0;    
    if (n == 0){
        for(int j = 1; j < 4; j++){
            for(int k = 0; k < 4; k++){
                result += pow(x_scaled, j - 1) * pow(y_scaled, k) * spline.f_coeffs[j][k] * float(j) / dx;
            }
        }
    } else if (n==1){
        for(int j = 0; j < 4; j++){
            for(int k = 1; k < 4; k++){
                result += pow(x_scaled, j) * pow(y_scaled, k - 1) * spline.f_coeffs[j][k] * float(k) / dy;
            }
        }
    } else {
        cout << "Interpolation derivative index must be 0 or 1 for a 2D function." << '\n';
        return 0.0;
    }
    
    return result;
}

double interp::eval_ddf(double x, double y, int n1, int n2, struct bicubic_spline & spline){
    int nx, ny;
    double dx, dy, x_scaled, y_scaled, result = 0.0;
    
    if(!same_cell(x, y, spline)){
        nx = find_segment(x, spline.x_vals, spline.length_x, spline.accel[0]);
        ny = find_segment(y, spline.y_vals, spline.length_y, spline.accel[1]);
        update_spline_coeffs(nx, ny, spline);
    } else {
        nx = spline.accel[0];
        ny = spline.accel[1];
    }
    
    dx = spline.x_vals[nx + 1] - spline.x_vals[nx];    x_scaled = (x - spline.x_vals[nx]) / dx;
    dy = spline.y_vals[ny + 1] - spline.y_vals[ny];    y_scaled = (y - spline.y_vals[ny]) / dy;

    result = 0.0;    
    if (n1 == 0 && n2 == 0){
        for(int j = 2; j < 4; j++){
            for(int k = 0; k < 4; k++){
                result += pow(x_scaled, j - 2) * pow(y_scaled, k) * spline.f_coeffs[j][k] * float(j * (j - 1)) / (dx * dx);
            }
        }
    } else if (n1 == 1 && n2 == 1){
        for(int j = 0; j < 4; j++){
            for(int k = 2; k < 4; k++){
                result += pow(x_scaled, j) * pow(y_scaled, k - 2) * spline.f_coeffs[j][k] * float(k * (k - 1)) / (dy * dy);
            }
        }
    } else if ((n1 == 0 && n2 == 1) || (n1 == 1 && n2 == 0)){
        for(int j = 1; j < 4; j++){
            for(int k = 1; k < 4; k++){
                result += pow(x_scaled, j - 1) * pow(y_scaled, k - 1) * spline.f_coeffs[j][k] * float(j * k) / (dx * dy);
            }
        }
    } else {
        cout << "Interpolation derivative index must be 0 or 1 for a 2D function." << '\n';
        return 0.0;
    }
    
    return result;
}

//---------------------------------------//
//-----------Three-Dimensional-----------//
//----------Hybrid Interpolation---------//
//---------------------------------------//
void interp::prep(struct hybrid_spline_3D & spline, int nx, int ny, int nz){
    spline.length_x = nx;   spline.accel[0] = 0;   spline.x_vals = new double [spline.length_x];
    spline.length_y = ny;   spline.accel[1] = 0;   spline.y_vals = new double [spline.length_y];
    spline.length_z = nz;   spline.accel[2] = 0;   spline.z_vals = new double [spline.length_z];
    
    spline.f_vals = new double** [spline.length_x];
    spline.f_slopes = new double** [spline.length_x];
    spline.dfdx_slopes = new double** [spline.length_x];
    spline.dfdy_slopes = new double** [spline.length_x];
    
    for(int nx = 0; nx < spline.length_x; nx++){
        spline.f_vals[nx] = new double* [spline.length_y];
        spline.f_slopes[nx] = new double* [spline.length_y];
        spline.dfdx_slopes[nx] = new double* [spline.length_y];
        spline.dfdy_slopes[nx] = new double* [spline.length_y];
        
        for(int ny = 0; ny < spline.length_y; ny++){
            spline.f_vals[nx][ny] = new double [spline.length_z];
            spline.f_slopes[nx][ny] = new double [spline.length_z];
            spline.dfdx_slopes[nx][ny] = new double [spline.length_z];
            spline.dfdy_slopes[nx][ny] = new double [spline.length_z];
        }
    }
}

void interp::set(struct hybrid_spline_3D & spline){
    int mx_up, mx_dn, my_up, my_dn;
    double dfdx [spline.length_z];
    double dfdy [spline.length_z];

    double dz1, dz2, ai, bi, ci, di;
    
    double new_c[spline.length_z - 1];
    double new_d[spline.length_z];
    
    for(int mx = 0; mx < spline.length_x; mx++){
        for(int my = 0; my < spline.length_y; my++){
            bi = 2.0 / (spline.z_vals[1] - spline.z_vals[0]);
            ci = 1.0 / (spline.z_vals[1] - spline.z_vals[0]);
            di = 3.0 * (spline.f_vals[mx][my][1] - spline.f_vals[mx][my][0]) / pow(spline.z_vals[1] - spline.z_vals[0], 2);
            
            new_c[0] = ci / bi;
            new_d[0] = di / bi;
            
            for(int i = 1; i < spline.length_z - 1; i++) {
                dz1 = spline.z_vals[i] - spline.z_vals[i - 1];
                dz2 = spline.z_vals[i + 1] - spline.z_vals[i];
            
                ai = 1.0 / dz1; bi = 2.0 * (1.0 / dz1 + 1.0 / dz2); ci = 1.0 / dz2;
                di = 3.0 * ((spline.f_vals[mx][my][i] - spline.f_vals[mx][my][i - 1]) / pow(dz1, 2)
                          + (spline.f_vals[mx][my][i + 1] - spline.f_vals[mx][my][i]) / pow(dz2, 2));

                new_c[i] = ci / (bi - new_c[i - 1] * ai);
                new_d[i] = (di - new_d[i - 1] * ai)/(bi - new_c[i - 1] * ai);
            }
            
            ai = 1.0 / (spline.z_vals[spline.length_z - 1] - spline.z_vals[spline.length_z - 2]);
            bi = 2.0 / (spline.z_vals[spline.length_z - 1] - spline.z_vals[spline.length_z - 2]);
            di = 3.0 * (spline.f_vals[mx][my][spline.length_z - 1] - spline.f_vals[mx][my][spline.length_z - 2])
                            / pow(spline.z_vals[spline.length_z - 1] - spline.z_vals[spline.length_z - 2], 2);
            
            new_d[spline.length_z-1] = (di - new_d[spline.length_z - 2] * ai) / (bi - new_c[spline.length_z - 2] * ai);
            
            spline.f_slopes[mx][my][spline.length_z - 1] = new_d[spline.length_z - 1];
            for(int i = spline.length_z - 2; i >= 0; i--){
                spline.f_slopes[mx][my][i] = new_d[i] - new_c[i] * spline.f_slopes[mx][my][i + 1];
            }
        }
    }
    
    for(int mx = 0; mx < spline.length_x; mx++){
        for(int my = 0; my < spline.length_y; my++){
            for(int mz = 0; mz < spline.length_z; mz++){
                mx_up = min(mx + 1, spline.length_x - 1);
                my_up = min(my + 1, spline.length_y - 1);
                
                mx_dn = max(mx - 1, 0);
                my_dn = max(my - 1, 0);
                
                dfdx[mz] = (spline.f_vals[mx_up][my][mz] - spline.f_vals[mx_dn][my][mz]) / (spline.x_vals[mx_up] - spline.x_vals[mx_dn]);
                dfdy[mz] = (spline.f_vals[mx][my_up][mz] - spline.f_vals[mx][my_dn][mz]) / (spline.y_vals[my_up] - spline.y_vals[my_dn]);
            }

            bi = 2.0 / (spline.z_vals[1] - spline.z_vals[0]);
            ci = 1.0 / (spline.z_vals[1] - spline.z_vals[0]);
            di = 3.0 * (dfdx[1] - dfdx[0]) / pow(spline.z_vals[1] - spline.z_vals[0], 2);
            
            new_c[0] = ci / bi;
            new_d[0] = di / bi;
   
            for(int i = 1; i < spline.length_z - 1; i++) {
                dz1 = spline.z_vals[i] - spline.z_vals[i - 1],
                dz2 = spline.z_vals[i + 1] - spline.z_vals[i];
            
                ai = 1.0 / dz1; bi = 2.0 * (1.0 / dz1 + 1.0 / dz2); ci = 1.0 / dz2;
                di = 3.0 * ((dfdx[i] - dfdx[i - 1]) / pow(dz1, 2)
                          + (dfdx[i + 1] - dfdx[i]) / pow(dz2, 2));
                
                new_c[i] = ci/(bi - new_c[i - 1] * ai);
                new_d[i] = (di - new_d[i-1] * ai) / (bi - new_c[i-1] * ai);
            }
            
            ai = 1.0 / (spline.z_vals[spline.length_z - 1] - spline.z_vals[spline.length_z - 2]);
            bi = 2.0 / (spline.z_vals[spline.length_z - 1] - spline.z_vals[spline.length_z - 2]);
            di = 3.0 * (dfdx[spline.length_z - 1] - dfdx[spline.length_z - 2])
                        / pow(spline.z_vals[spline.length_z-1] - spline.z_vals[spline.length_z - 2], 2);
            
            new_d[spline.length_z - 1] = (di - new_d[spline.length_z - 2] * ai) / (bi - new_c[spline.length_z - 2] * ai);
            spline.dfdx_slopes[mx][my][spline.length_z - 1] = new_d[spline.length_z - 1];
            for(int i = spline.length_z - 2; i >= 0; i--){
                spline.dfdx_slopes[mx][my][i] = new_d[i] - new_c[i] * spline.dfdx_slopes[mx][my][i + 1];
            }
        }
    }
    
    for(int mx = 0; mx < spline.length_x; mx++){
        for(int my = 0; my < spline.length_y; my++){
            for(int mz = 0; mz < spline.length_z; mz++){
                mx_up = min(mx + 1, spline.length_x - 1);
                my_up = min(my + 1, spline.length_y - 1);
                
                mx_dn = max(mx - 1, 0);
                my_dn = max(my - 1, 0);
                
                dfdx[mz] = (spline.f_vals[mx_up][my][mz] - spline.f_vals[mx_dn][my][mz]) / (spline.x_vals[mx_up] - spline.x_vals[mx_dn]);
                dfdy[mz] = (spline.f_vals[mx][my_up][mz] - spline.f_vals[mx][my_dn][mz]) / (spline.y_vals[my_up] - spline.y_vals[my_dn]);
            }

            bi = 2.0 / (spline.z_vals[1] - spline.z_vals[0]);
            ci = 1.0 / (spline.z_vals[1] - spline.z_vals[0]);
            di = 3.0 * (dfdy[1] - dfdy[0]) / pow(spline.z_vals[1] - spline.z_vals[0], 2);
            
            new_c[0] = ci / bi;
            new_d[0] = di / bi;

            for(int i = 1; i < spline.length_z - 1; i++) {
                dz1 = spline.z_vals[i] - spline.z_vals[i - 1];
                dz2 = spline.z_vals[i + 1] - spline.z_vals[i];
            
                ai = 1.0 / dz1; bi = 2.0 * (1.0 / dz1 + 1.0 / dz2); ci = 1.0 / dz2;
                di = 3.0 * ((dfdy[i] - dfdy[i - 1]) / pow(dz1, 2)
                          + (dfdy[i + 1] - dfdy[i]) / pow(dz2, 2));
                
                new_c[i] = ci / (bi - new_c[i - 1] * ai);
                new_d[i] = (di - new_d[i - 1] * ai) / (bi - new_c[i - 1] * ai);
            }
            
            ai = 1.0/(spline.z_vals[spline.length_z - 1] - spline.z_vals[spline.length_z - 2]);
            bi = 2.0/(spline.z_vals[spline.length_z - 1] - spline.z_vals[spline.length_z - 2]);
            di = 3.0 * (dfdy[spline.length_z - 1] - dfdy[spline.length_z - 2])
                    / pow(spline.z_vals[spline.length_z - 1] - spline.z_vals[spline.length_z - 2], 2);
            
            new_d[spline.length_z - 1] = (di - new_d[spline.length_z - 2] * ai) / (bi - new_c[spline.length_z - 2] * ai);
            spline.dfdy_slopes[mx][my][spline.length_z - 1] = new_d[spline.length_z - 1];
            for(int i = spline.length_z - 2; i >= 0; i--){
                spline.dfdy_slopes[mx][my][i] = new_d[i] - new_c[i] * spline.dfdy_slopes[mx][my][i + 1];
            }
        }
    }
}

void interp::clear(struct hybrid_spline_3D & spline){
    delete spline.x_vals;
    delete spline.y_vals;
    delete spline.z_vals;
    
    for(int nx = 0; nx < spline.length_x; nx++){
        for(int ny = 0; ny < spline.length_y; ny++){
            delete spline.f_vals[nx][ny];
            delete spline.f_slopes[nx][ny];
            delete spline.dfdx_slopes[nx][ny];
            delete spline.dfdy_slopes[nx][ny];
        }
        delete spline.f_vals[nx];
        delete spline.f_slopes[nx];
        delete spline.dfdx_slopes[nx];
        delete spline.dfdy_slopes[nx];
    }
    delete spline.f_vals;
    delete spline.f_slopes;
    delete spline.dfdx_slopes;
    delete spline.dfdy_slopes;
}

// Evaluate the functions on the nodes
double interp::eval_node_f(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    double dz, X, df, A, B;

    dz = spline.z_vals[kz + 1] - spline.z_vals[kz];
    X = (z - spline.z_vals[kz]) / dz;
    
    df = spline.f_vals[kx][ky][kz + 1] - spline.f_vals[kx][ky][kz];
    A = spline.f_slopes[kx][ky][kz] * dz - df;
    B = -spline.f_slopes[kx][ky][kz + 1] * dz + df;
    
    return spline.f_vals[kx][ky][kz] + X * (df + (1.0 - X) * (A * (1.0 - X) + B * X));
}

double interp::eval_node_dfdx(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    int kx_up, kx_dn;
    double dfdx_kz, dfdx_kzp1, df, dz, X, A, B;
    
    kx_up = min(kx + 1, spline.length_x - 1);
    kx_dn = max(kx - 1, 0);
    
    dfdx_kz = (spline.f_vals[kx_up][ky][kz] - spline.f_vals[kx_dn][ky][kz])/(spline.x_vals[kx_up] - spline.x_vals[kx_dn]);
    dfdx_kzp1 = (spline.f_vals[kx_up][ky][kz + 1] - spline.f_vals[kx_dn][ky][kz + 1])/(spline.x_vals[kx_up] - spline.x_vals[kx_dn]);
    
    dz = spline.z_vals[kz + 1] - spline.z_vals[kz];
    X = (z - spline.z_vals[kz]) / dz; 
    
    df = dfdx_kzp1 - dfdx_kz; 
    A = spline.dfdx_slopes[kx][ky][kz] * dz - df;
    B = -spline.dfdx_slopes[kx][ky][kz + 1] * dz + df;
    
    return dfdx_kz + X * (df + (1.0 - X) * (A * (1.0 - X) + B * X));
}

double interp::eval_node_dfdy(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    int ky_up, ky_dn;
    double dfdy_kz, dfdy_kzp1, df, dz, X, A, B;
    
    ky_up = min(ky + 1, spline.length_y - 1);
    ky_dn = max(ky - 1, 0);
    
    dfdy_kz =   (spline.f_vals[kx][ky_up][kz] - spline.f_vals[kx][ky_dn][kz]) / (spline.y_vals[ky_up] - spline.y_vals[ky_dn]);
    dfdy_kzp1 = (spline.f_vals[kx][ky_up][kz + 1] - spline.f_vals[kx][ky_dn][kz + 1]) / (spline.y_vals[ky_up] - spline.y_vals[ky_dn]);
    
    dz = spline.z_vals[kz + 1] - spline.z_vals[kz];
    X = (z - spline.z_vals[kz]) / dz; 
    
    df = dfdy_kzp1 - dfdy_kz;
    A = spline.dfdx_slopes[kx][ky][kz] * dz - df;
    B = -spline.dfdx_slopes[kx][ky][kz + 1] * dz + df;
    
    return dfdy_kz + X * (df + (1.0 - X) * (A * (1.0 - X) + B * X));
}

double interp::eval_node_dfdz(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    double dz, X, df, A, B;

    dz = spline.z_vals[kz + 1] - spline.z_vals[kz];
    X = (z - spline.z_vals[kz]) / dz;
    
    df = spline.f_vals[kx][ky][kz + 1] - spline.f_vals[kx][ky][kz];
    A = spline.f_slopes[kx][ky][kz] * dz - df;
    B = -spline.f_slopes[kx][ky][kz + 1] * dz + df;
    
    return df / dz + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X) / dz + X * (1.0 - X) * (B - A) / dz;
}

double interp::eval_node_ddfdxdz(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    int kx_up, kx_dn;
    double dfdx_kz, dfdx_kzp1, df, dz, X, A, B;
    
    kx_up = min(kx + 1, spline.length_x - 1);
    kx_dn = max(kx - 1, 0);
    
    dfdx_kz = (spline.f_vals[kx_up][ky][kz] - spline.f_vals[kx_dn][ky][kz])/(spline.x_vals[kx_up] - spline.x_vals[kx_dn]);
    dfdx_kzp1 = (spline.f_vals[kx_up][ky][kz + 1] - spline.f_vals[kx_dn][ky][kz + 1])/(spline.x_vals[kx_up] - spline.x_vals[kx_dn]);
    
    dz = spline.z_vals[kz + 1] - spline.z_vals[kz];
    X = (z - spline.z_vals[kz]) / dz;
    
    df = dfdx_kzp1 - dfdx_kz;
    A = spline.dfdx_slopes[kx][ky][kz] * dz - df;
    B = -spline.dfdx_slopes[kx][ky][kz + 1] * dz + df;
    
    return df/dz + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X) / dz + X * (1.0 - X) * (B - A)/ dz;
}

double interp::eval_node_ddfdydz(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    int ky_up, ky_dn;
    double dfdy_kz, dfdy_kzp1, df, dz, X, A, B;
    
    ky_up = min(ky + 1, spline.length_y - 1);
    ky_dn = max(ky - 1, 0);
    
    dfdy_kz =   (spline.f_vals[kx][ky_up][kz] - spline.f_vals[kx][ky_dn][kz]) / (spline.y_vals[ky_up] - spline.y_vals[ky_dn]);
    dfdy_kzp1 = (spline.f_vals[kx][ky_up][kz + 1] - spline.f_vals[kx][ky_dn][kz + 1]) / (spline.y_vals[ky_up] - spline.y_vals[ky_dn]);
    
    dz = spline.z_vals[kz + 1] - spline.z_vals[kz];
    X = (z - spline.z_vals[kz]) / dz;

    df = dfdy_kzp1 - dfdy_kz;
    A = spline.dfdx_slopes[kx][ky][kz] * dz - df;
    B = -spline.dfdx_slopes[kx][ky][kz + 1] * dz + df;
    
    return df / dz + (1.0 - 2.0 * X) * (A * (1.0 - X) + B * X) / dz + X * (1.0 - X) * (B - A) / dz;
}

double interp::eval_node_ddfdzdz(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    double dz, X, df, A, B;

    dz = spline.z_vals[kz + 1] - spline.z_vals[kz];
    X = (z - spline.z_vals[kz]) / dz;
    
    df = spline.f_vals[kx][ky][kz + 1] - spline.f_vals[kx][ky][kz];
    A = spline.f_slopes[kx][ky][kz] * dz - df;
    B = -spline.f_slopes[kx][ky][kz + 1] * dz + df;
    
    return 2.0 * (B - 2.0 * A + (A - B) * 3.0 * X) / pow(dz, 2);
}

// Finite difference derivatives for defining the bicubic spline corners
double interp::finite_diff_dfdx(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    int kx_up, kx_dn;

    kx_up = min(kx + 1, spline.length_x - 1);
    kx_dn = max(kx - 1, 0);
    
    return (eval_node_f(z, spline, kx_up, ky, kz) - eval_node_f(z, spline, kx_dn, ky, kz)) / (spline.x_vals[kx_up] - spline.x_vals[kx_dn]);
}

double interp::finite_diff_dfdy(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    int ky_up, ky_dn;
    
    ky_up = min(ky + 1, spline.length_y - 1);
    ky_dn = max(ky - 1, 0);
    
    return (eval_node_f(z, spline, kx, ky_up, kz) - eval_node_f(z, spline, kx, ky_dn, kz)) / (spline.y_vals[ky_up] - spline.y_vals[ky_dn]);
}

double interp::finite_diff_ddfdxdx(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    int kx_up, kx_dn;
    
    kx_up = min(kx + 1, spline.length_x - 1);
    kx_dn = max(kx - 1, 0);
    
    return (eval_node_dfdx(z, spline, kx_up, ky, kz) - eval_node_dfdx(z, spline, kx_dn, ky, kz)) / (spline.x_vals[kx_up] - spline.x_vals[kx_dn]);
}

double interp::finite_diff_ddfdydy(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    int ky_up, ky_dn;
    
    ky_up = min(ky + 1, spline.length_y - 1);
    ky_dn = max(ky - 1, 0);
    
    return (eval_node_dfdy(z, spline, kx, ky_up, kz) - eval_node_dfdy(z, spline, kx, ky_dn, kz)) / (spline.y_vals[ky_up] - spline.y_vals[ky_dn]);
}

double interp::finite_diff_ddfdxdy(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    int kx_up, ky_up, kx_dn, ky_dn;

    kx_up = min(kx + 1, spline.length_x - 1);
    ky_up = min(ky + 1, spline.length_y - 1);
    
    kx_dn = max(kx - 1, 0);
    ky_dn = max(ky - 1, 0);

    return (eval_node_f(z, spline, kx_up, ky_up, kz) - eval_node_f(z, spline, kx_up, ky_dn, kz) - eval_node_f(z, spline, kx_dn, ky_up, kz) + eval_node_f(z, spline, kx_dn, ky_dn, kz))
                    /((spline.x_vals[kx_up] - spline.x_vals[kx_dn]) * (spline.y_vals[ky_up] - spline.y_vals[ky_dn]));
}

double interp::finite_diff_ddfdxdz(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    int kx_up, kx_dn;

    kx_up = min(kx + 1, spline.length_x - 1);
    kx_dn = max(kx - 1, 0);
    
    return (eval_node_dfdz(z, spline, kx_up, ky, kz) - eval_node_dfdz(z, spline, kx_dn, ky, kz)) / (spline.x_vals[kx_up] - spline.x_vals[kx_dn]);
}

double interp::finite_diff_ddfdydz(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    int ky_up, ky_dn;
    
    ky_up = min(ky + 1, spline.length_y - 1);
    ky_dn = max(ky - 1, 0);
    
    return (eval_node_dfdz(z, spline, kx, ky_up, kz) - eval_node_dfdz(z, spline, kx, ky_dn, kz)) / (spline.y_vals[ky_up] - spline.y_vals[ky_dn]);
}


double interp::finite_diff_dddfdxdxdy(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    int kx_up, ky_up, kx_dn, ky_dn;

    kx_up = min(kx + 1, spline.length_x - 1);
    ky_up = min(ky + 1, spline.length_y - 1);
    
    kx_dn = max(kx - 1, 0);
    ky_dn = max(ky - 1, 0);
    
    return (eval_node_dfdx(z, spline, kx_up, ky_up, kz) - eval_node_dfdx(z, spline, kx_up, ky_dn, kz)
            - eval_node_dfdx(z, spline, kx_dn, ky_up, kz) + eval_node_dfdx(z, spline, kx_dn, ky_dn, kz))
                /((spline.x_vals[kx_up] - spline.x_vals[kx_dn]) * (spline.y_vals[ky_up] - spline.y_vals[ky_dn]));
}

double interp::finite_diff_dddfdxdydy(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    int kx_up, ky_up, kx_dn, ky_dn;

    kx_up = min(kx + 1, spline.length_x - 1);
    ky_up = min(ky + 1, spline.length_y - 1);
    
    kx_dn = max(kx - 1, 0);
    ky_dn = max(ky - 1, 0);

    return (eval_node_dfdy(z, spline, kx_up, ky_up, kz) - eval_node_dfdy(z, spline, kx_up, ky_dn, kz)
            - eval_node_dfdy(z, spline, kx_dn, ky_up, kz) + eval_node_dfdy(z, spline, kx_dn, ky_dn, kz))
                /((spline.x_vals[kx_up] - spline.x_vals[kx_dn])*(spline.y_vals[ky_up] - spline.y_vals[ky_dn]));
}

double interp::finite_diff_dddfdxdydz(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    int kx_up, ky_up, kx_dn, ky_dn;

    kx_up = min(kx + 1, spline.length_x - 1);
    ky_up = min(ky + 1, spline.length_y - 1);
    
    kx_dn = max(kx - 1, 0);
    ky_dn = max(ky - 1, 0);
    
    return (eval_node_dfdz(z, spline, kx_up, ky_up, kz) - eval_node_dfdz(z, spline, kx_up, ky_dn, kz)
            - eval_node_dfdz(z, spline, kx_dn, ky_up, kz) + eval_node_dfdz(z, spline, kx_dn, ky_dn, kz))
                /((spline.x_vals[kx_up] - spline.x_vals[kx_dn]) * (spline.y_vals[ky_up] - spline.y_vals[ky_dn]));
}

double interp::finite_diff_dddfdxdzdz(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    int kx_up, kx_dn;

    kx_up = min(kx + 1, spline.length_x - 1);
    kx_dn = max(kx - 1, 0);

    return (eval_node_ddfdzdz(z, spline, kx_up, ky, kz) - eval_node_ddfdzdz(z, spline, kx_dn, ky, kz)) / (spline.x_vals[kx_up] - spline.x_vals[kx_dn]);
}

double interp::finite_diff_dddfdydzdz(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    int ky_up, ky_dn;
    
    ky_up = min(ky + 1, spline.length_y - 1);
    ky_dn = max(ky - 1, 0);
    
    return (eval_node_ddfdzdz(z, spline, kx, ky_up, kz) - eval_node_ddfdzdz(z, spline, kx, ky_dn, kz)) / (spline.y_vals[ky_up] - spline.y_vals[ky_dn]);

}

double interp::finite_diff_ddddfdxdydzdz(double z, struct hybrid_spline_3D spline, int kx, int ky, int kz){
    int kx_up, ky_up, kx_dn, ky_dn;

    kx_up = min(kx + 1, spline.length_x - 1);
    ky_up = min(ky + 1, spline.length_y - 1);
    
    kx_dn = max(kx - 1, 0);
    ky_dn = max(ky - 1, 0);
        
    return (eval_node_ddfdzdz(z, spline, kx_up, ky_up, kz) - eval_node_ddfdzdz(z, spline, kx_up, ky_dn, kz)
            - eval_node_ddfdzdz(z, spline, kx_dn, ky_up, kz) + eval_node_ddfdzdz(z, spline, kx_dn, ky_dn, kz))
                /((spline.x_vals[kx_up] - spline.x_vals[kx_dn]) * (spline.y_vals[ky_up] - spline.y_vals[ky_dn]));
}

double interp::eval_f(double x, double y, double z, struct hybrid_spline_3D & spline){
    int kx, ky, kz;
    double dx, dy, x_scaled, y_scaled, X_vec[16], A_vec[16], result;
    
    kx = find_segment(x, spline.x_vals, spline.length_x, spline.accel[0]);
    ky = find_segment(y, spline.y_vals, spline.length_y, spline.accel[1]);
    kz = find_segment(z, spline.z_vals, spline.length_z, spline.accel[2]);

    dx = spline.x_vals[kx + 1] - spline.x_vals[kx];  x_scaled = (x - spline.x_vals[kx]) / dx;
    dy = spline.y_vals[ky + 1] - spline.y_vals[ky];  y_scaled = (y - spline.y_vals[ky]) / dy;
    
    X_vec[0] = eval_node_f(z, spline, kx, ky, kz);                  X_vec[8] =  finite_diff_dfdy(z, spline, kx, ky, kz) * dy;
    X_vec[1] = eval_node_f(z, spline, kx + 1, ky, kz);              X_vec[9] =  finite_diff_dfdy(z, spline, kx + 1, ky, kz) * dy;
    X_vec[2] = eval_node_f(z, spline, kx, ky + 1, kz);              X_vec[10] = finite_diff_dfdy(z, spline, kx, ky + 1, kz) * dy;
    X_vec[3] = eval_node_f(z, spline, kx + 1, ky + 1, kz);          X_vec[11] = finite_diff_dfdy(z, spline, kx + 1, ky + 1, kz) * dy;
    
    X_vec[4] = finite_diff_dfdx(z, spline, kx, ky, kz) * dx;        X_vec[12] = finite_diff_ddfdxdy(z, spline, kx, ky, kz) * dx * dy;
    X_vec[5] = finite_diff_dfdx(z, spline, kx + 1, ky, kz) * dx;    X_vec[13] = finite_diff_ddfdxdy(z, spline, kx + 1, ky, kz) * dx * dy;
    X_vec[6] = finite_diff_dfdx(z, spline, kx, ky + 1, kz) * dx;    X_vec[14] = finite_diff_ddfdxdy(z, spline, kx, ky + 1, kz) * dx * dy;
    X_vec[7] = finite_diff_dfdx(z, spline, kx + 1, ky+1, kz) * dx;  X_vec[15] = finite_diff_ddfdxdy(z, spline, kx + 1, ky + 1, kz) * dx * dy;
    
    for(int j = 0; j < 16; j++){
        A_vec[j] = 0.0;
        for(int k = 0; k < 16; k++){
            A_vec[j] += bicubic_conversion_matrix[j][k] * X_vec[k];
        }
    }
	
    result = 0.0;
    for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            result += A_vec[k2 * 4 + k1] * pow(x_scaled, k1) * pow(y_scaled, k2);
        }
    }
	
    return result;
}

double interp::eval_df(double x, double y, double z, int n, struct hybrid_spline_3D & spline){
    int kx, ky, kz;
    double dx, dy, x_scaled, y_scaled, X_vec[16], A_vec[16], result;
    
    kx = find_segment(x, spline.x_vals, spline.length_x, spline.accel[0]);
    ky = find_segment(y, spline.y_vals, spline.length_y, spline.accel[1]);
    kz = find_segment(z, spline.z_vals, spline.length_z, spline.accel[2]);
    
    dx = spline.x_vals[kx + 1] - spline.x_vals[kx];  x_scaled = (x - spline.x_vals[kx]) / dx;
    dy = spline.y_vals[ky + 1] - spline.y_vals[ky];  y_scaled = (y - spline.y_vals[ky]) / dy;
    
    if(n == 0 || n == 1){
        // df/dx or df/dy from d/dx or d/dy of bicubic interpolation of f
        X_vec[0] = eval_node_f(z, spline, kx, ky, kz);                      X_vec[8] =  finite_diff_dfdy(z, spline, kx, ky, kz) * dy;
        X_vec[1] = eval_node_f(z, spline, kx + 1, ky, kz);                  X_vec[9] =  finite_diff_dfdy(z, spline, kx + 1, ky, kz) * dy;
        X_vec[2] = eval_node_f(z, spline, kx, ky + 1, kz);                  X_vec[10] = finite_diff_dfdy(z, spline, kx, ky + 1, kz) * dy;
        X_vec[3] = eval_node_f(z, spline, kx + 1, ky + 1, kz);              X_vec[11] = finite_diff_dfdy(z, spline, kx + 1, ky + 1, kz) * dy;
        
        X_vec[4] = finite_diff_dfdx(z, spline, kx, ky, kz) * dx;            X_vec[12] = finite_diff_ddfdxdy(z, spline, kx, ky, kz) * dx * dy;
        X_vec[5] = finite_diff_dfdx(z, spline, kx + 1, ky, kz) * dx;        X_vec[13] = finite_diff_ddfdxdy(z, spline, kx + 1, ky, kz) * dx * dy;
        X_vec[6] = finite_diff_dfdx(z, spline, kx, ky + 1, kz) * dx;        X_vec[14] = finite_diff_ddfdxdy(z, spline, kx, ky + 1, kz) * dx * dy;
        X_vec[7] = finite_diff_dfdx(z, spline, kx + 1, ky+1, kz) * dx;      X_vec[15] = finite_diff_ddfdxdy(z, spline, kx + 1, ky + 1, kz) * dx * dy;
    } else {
        // df/dz from bicubic interpolation of df/dz
        X_vec[0] = eval_node_dfdz(z, spline, kx, ky, kz);                   X_vec[8] =  finite_diff_ddfdydz(z, spline, kx, ky, kz) * dy;
        X_vec[1] = eval_node_dfdz(z, spline, kx + 1, ky, kz);               X_vec[9] =  finite_diff_ddfdydz(z, spline, kx + 1, ky, kz) * dy;
        X_vec[2] = eval_node_dfdz(z, spline, kx, ky + 1, kz);               X_vec[10] = finite_diff_ddfdydz(z, spline, kx, ky + 1, kz) * dy;
        X_vec[3] = eval_node_dfdz(z, spline, kx + 1, ky + 1, kz);           X_vec[11] = finite_diff_ddfdydz(z, spline, kx + 1, ky + 1, kz) * dy;
        
        X_vec[4] = finite_diff_ddfdxdz(z, spline, kx, ky, kz) * dx;         X_vec[12] = finite_diff_dddfdxdydz(z, spline, kx, ky, kz) * dx * dy;
        X_vec[5] = finite_diff_ddfdxdz(z, spline, kx + 1, ky ,kz) * dx;     X_vec[13] = finite_diff_dddfdxdydz(z, spline, kx + 1, ky, kz) * dx * dy;
        X_vec[6] = finite_diff_ddfdxdz(z, spline, kx, ky + 1, kz) * dx;     X_vec[14] = finite_diff_dddfdxdydz(z, spline, kx, ky + 1, kz) * dx * dy;
        X_vec[7] = finite_diff_ddfdxdz(z, spline, kx + 1, ky + 1, kz) * dx; X_vec[15] = finite_diff_dddfdxdydz(z, spline, kx + 1, ky + 1, kz) * dx * dy;
    }
    
    for(int j = 0; j < 16; j++){
        A_vec[j] = 0;
        for(int k = 0; k < 16; k++){
            A_vec[j] += bicubic_conversion_matrix[j][k] * X_vec[k];
        }
    }

    result = 0.0;
    if(n == 0){
        for(int k1 = 1; k1 < 4; k1++){
            for(int k2 = 0; k2 < 4; k2++){
                result += A_vec[k2 * 4 + k1] * pow(x_scaled, k1 - 1) * pow(y_scaled, k2) * float(k1) / dx;
            }
        }
    } else if(n == 1){
        for(int k1 = 0; k1 < 4; k1++){
            for(int k2 = 1; k2 < 4; k2++){
                result += A_vec[k2 * 4 + k1] * pow(x_scaled, k1) * pow(y_scaled, k2 - 1) * float(k2) / dy;
            }
        }
    } else {
        for(int k1 = 0; k1 < 4; k1++){
            for(int k2 = 0; k2 < 4; k2++){
                result += A_vec[k2 * 4 + k1] * pow(x_scaled, k1) * pow(y_scaled, k2);
            }
        }
    }
    
    return result;
}

double interp::eval_ddf(double x, double y, double z, int n1, int n2, struct hybrid_spline_3D & spline){
    int kx, ky, kz;
    double dx, dy, x_scaled, y_scaled, X_vec[16], A_vec[16], result;
    
    kx = find_segment(x, spline.x_vals, spline.length_x, spline.accel[0]);
    ky = find_segment(y, spline.y_vals, spline.length_y, spline.accel[1]);
    kz = find_segment(z, spline.z_vals, spline.length_z, spline.accel[2]);
    
    dx = spline.x_vals[kx + 1] - spline.x_vals[kx];  x_scaled = (x - spline.x_vals[kx]) / dx;
    dy = spline.y_vals[ky + 1] - spline.y_vals[ky];  y_scaled = (y - spline.y_vals[ky]) / dy;
    
    if(n1 == 0 && n2 == 0){
        // ddf/dxdx from d/dx of bicubic interpolation of df/dx
        X_vec[0] = eval_node_dfdx(z, spline, kx, ky, kz);                       X_vec[8] =  finite_diff_ddfdxdy(z, spline, kx, ky, kz) * dy;
        X_vec[1] = eval_node_dfdx(z, spline, kx + 1, ky, kz);                   X_vec[9] =  finite_diff_ddfdxdy(z, spline, kx + 1, ky, kz) * dy;
        X_vec[2] = eval_node_dfdx(z, spline, kx, ky + 1, kz);                   X_vec[10] = finite_diff_ddfdxdy(z, spline, kx, ky + 1, kz) * dy;
        X_vec[3] = eval_node_dfdx(z, spline, kx + 1, ky + 1, kz);               X_vec[11] = finite_diff_ddfdxdy(z, spline, kx + 1, ky + 1, kz) * dy;
        
        X_vec[4] = finite_diff_ddfdxdx(z, spline, kx, ky, kz) * dx;             X_vec[12] = finite_diff_dddfdxdxdy(z, spline, kx, ky, kz) * dx * dy;
        X_vec[5] = finite_diff_ddfdxdx(z, spline, kx + 1, ky, kz) * dx;         X_vec[13] = finite_diff_dddfdxdxdy(z, spline, kx + 1, ky, kz) * dx * dy;
        X_vec[6] = finite_diff_ddfdxdx(z, spline, kx, ky + 1, kz) * dx;         X_vec[14] = finite_diff_dddfdxdxdy(z, spline, kx, ky + 1, kz) * dx * dy;
        X_vec[7] = finite_diff_ddfdxdx(z, spline, kx + 1, ky + 1, kz) * dx;     X_vec[15] = finite_diff_dddfdxdxdy(z, spline, kx + 1, ky + 1, kz) * dx * dy;
    } else if(n1 == 1 && n2 == 1){
        // ddf/dydy from d/dy of bicubic interpolation of df/dy
        X_vec[0] = eval_node_dfdy(z, spline, kx, ky, kz);                       X_vec[8] =  finite_diff_ddfdydy(z, spline, kx, ky,kz) * dy;
        X_vec[1] = eval_node_dfdy(z, spline, kx + 1, ky, kz);                   X_vec[9] =  finite_diff_ddfdydy(z, spline, kx + 1, ky,kz) * dy;
        X_vec[2] = eval_node_dfdy(z, spline, kx, ky + 1, kz);                   X_vec[10] = finite_diff_ddfdydy(z, spline, kx, ky + 1,kz) * dy;
        X_vec[3] = eval_node_dfdy(z, spline, kx + 1, ky + 1, kz);               X_vec[11] = finite_diff_ddfdydy(z, spline, kx + 1, ky + 1,kz) * dy;
        
        X_vec[4] = finite_diff_ddfdxdy(z, spline, kx, ky, kz) * dx;             X_vec[12] = finite_diff_dddfdxdydy(z, spline, kx, ky, kz) * dx * dy;
        X_vec[5] = finite_diff_ddfdxdy(z, spline, kx + 1, ky, kz) * dx;         X_vec[13] = finite_diff_dddfdxdydy(z, spline, kx + 1, ky, kz) * dx * dy;
        X_vec[6] = finite_diff_ddfdxdy(z, spline, kx, ky + 1, kz) * dx;         X_vec[14] = finite_diff_dddfdxdydy(z, spline, kx, ky + 1, kz) * dx * dy;
        X_vec[7] = finite_diff_ddfdxdy(z, spline, kx + 1, ky + 1, kz) * dx;     X_vec[15] = finite_diff_dddfdxdydy(z, spline, kx + 1, ky + 1, kz) * dx * dy;
    } else if(n1 == 2 && n2 == 2){
        // ddf/dzdz from bicubic interpolation of ddf/dzdz
        X_vec[0] = eval_node_ddfdzdz(z, spline, kx, ky, kz);                    X_vec[8] =  finite_diff_dddfdydzdz(z, spline, kx, ky, kz) * dy;
        X_vec[1] = eval_node_ddfdzdz(z, spline, kx + 1, ky, kz);                X_vec[9] =  finite_diff_dddfdydzdz(z, spline, kx + 1, ky, kz) * dy;
        X_vec[2] = eval_node_ddfdzdz(z, spline, kx, ky + 1, kz);                X_vec[10] = finite_diff_dddfdydzdz(z, spline, kx, ky + 1, kz) * dy;
        X_vec[3] = eval_node_ddfdzdz(z, spline, kx + 1, ky + 1, kz);            X_vec[11] = finite_diff_dddfdydzdz(z, spline, kx + 1, ky + 1, kz) * dy;
        
        X_vec[4] = finite_diff_dddfdxdzdz(z, spline, kx, ky, kz) * dx;          X_vec[12] = finite_diff_ddddfdxdydzdz(z, spline, kx, ky, kz) * dx * dy;
        X_vec[5] = finite_diff_dddfdxdzdz(z, spline, kx + 1, ky, kz) * dx;      X_vec[13] = finite_diff_ddddfdxdydzdz(z, spline, kx + 1, ky, kz) * dx * dy;
        X_vec[6] = finite_diff_dddfdxdzdz(z, spline, kx, ky + 1, kz) * dx;      X_vec[14] = finite_diff_ddddfdxdydzdz(z, spline, kx, ky + 1, kz) * dx * dy;
        X_vec[7] = finite_diff_dddfdxdzdz(z, spline, kx + 1, ky + 1, kz) * dx;  X_vec[15] = finite_diff_ddddfdxdydzdz(z, spline, kx + 1, ky + 1, kz) * dx * dy;
    } else if((n1 == 0 && n2 == 1) || (n1 == 1 && n2 == 0)){
        // ddf/dxdy from dd/dxdy of bicubic interpolation of f
        X_vec[0] = eval_node_f(z, spline, kx, ky, kz);                          X_vec[8] =  finite_diff_dfdy(z, spline, kx, ky, kz) * dy;
        X_vec[1] = eval_node_f(z, spline, kx + 1, ky, kz);                      X_vec[9] =  finite_diff_dfdy(z, spline, kx + 1, ky, kz) * dy;
        X_vec[2] = eval_node_f(z, spline, kx, ky + 1, kz);                      X_vec[10] = finite_diff_dfdy(z, spline, kx, ky + 1, kz) * dy;
        X_vec[3] = eval_node_f(z, spline, kx + 1, ky + 1, kz);                  X_vec[11] = finite_diff_dfdy(z, spline, kx + 1, ky + 1, kz) * dy;
        
        X_vec[4] = finite_diff_dfdx(z, spline, kx, ky, kz) * dx;                X_vec[12] = finite_diff_ddfdxdy(z, spline, kx, ky, kz) * dx * dy;
        X_vec[5] = finite_diff_dfdx(z, spline, kx + 1, ky, kz) * dx;            X_vec[13] = finite_diff_ddfdxdy(z, spline, kx + 1, ky, kz) * dx * dy;
        X_vec[6] = finite_diff_dfdx(z, spline, kx, ky + 1, kz) * dx;            X_vec[14] = finite_diff_ddfdxdy(z, spline, kx, ky + 1, kz) * dx * dy;
        X_vec[7] = finite_diff_dfdx(z, spline, kx + 1, ky+1, kz) * dx;          X_vec[15] = finite_diff_ddfdxdy(z, spline, kx + 1, ky + 1, kz) * dx * dy;
        
    } else if((n1 == 0 && n2 == 2) || (n1 == 2 && n2 == 0) || (n1 == 1 && n2 == 2) || (n1 == 2 && n2 == 1)){
        // ddf/dxdz or ddf/dydz from d/dx or d/dy of bicubic interpolation of df/dz
        X_vec[0] = eval_node_dfdz(z, spline, kx, ky, kz);                       X_vec[8] =  finite_diff_ddfdydz(z, spline, kx, ky, kz) * dy;
        X_vec[1] = eval_node_dfdz(z, spline, kx + 1, ky, kz);                   X_vec[9] =  finite_diff_ddfdydz(z, spline, kx + 1, ky, kz) * dy;
        X_vec[2] = eval_node_dfdz(z, spline, kx, ky + 1, kz);                   X_vec[10] = finite_diff_ddfdydz(z, spline, kx, ky + 1, kz) * dy;
        X_vec[3] = eval_node_dfdz(z, spline, kx + 1, ky + 1, kz);               X_vec[11] = finite_diff_ddfdydz(z, spline, kx + 1, ky + 1, kz) * dy;
        
        X_vec[4] = finite_diff_ddfdxdz(z, spline, kx, ky, kz) * dx;             X_vec[12] = finite_diff_dddfdxdydz(z, spline, kx, ky, kz) * dx * dy;
        X_vec[5] = finite_diff_ddfdxdz(z, spline, kx + 1, ky ,kz) * dx;         X_vec[13] = finite_diff_dddfdxdydz(z, spline, kx + 1, ky, kz) * dx * dy;
        X_vec[6] = finite_diff_ddfdxdz(z, spline, kx, ky + 1, kz) * dx;         X_vec[14] = finite_diff_dddfdxdydz(z, spline, kx, ky + 1, kz) * dx * dy;
        X_vec[7] = finite_diff_ddfdxdz(z, spline, kx + 1, ky + 1, kz) * dx;     X_vec[15] = finite_diff_dddfdxdydz(z, spline, kx + 1, ky + 1, kz) * dx * dy;
    }
    
    for(int j = 0; j < 16; j++){
        A_vec[j] = 0.0;
        for(int k = 0; k < 16; k++){
            A_vec[j] += bicubic_conversion_matrix[j][k] * X_vec[k];
        }
    }

    result = 0.0;
    if (n1 == 2 && n2 == 2) {
        for(int k1 = 0; k1 < 4; k1++){
            for(int k2 = 0; k2 < 4; k2++){
                result += A_vec[k2 * 4 + k1] * pow(x_scaled, k1) * pow(y_scaled, k2);
            }
        }
    } else if((n1 == 0 && n2 == 0) || (n1 == 0 && n2 == 2) || (n1 == 2 && n2 == 0)){
        for(int k1 = 1; k1 < 4; k1++){
            for(int k2 = 0; k2 < 4; k2++){
                result += A_vec[k2 * 4 + k1] * pow(x_scaled, k1 - 1) * pow(y_scaled, k2) * float(k1) / dx;
            }
        }
    } else if((n1 == 1 && n2 == 1) || (n1 == 1 && n2 == 2) || (n1 == 2 && n2 == 1)){
        for(int k1 = 0; k1 < 4; k1++){
            for(int k2 = 1; k2 < 4; k2++){
                result += A_vec[k2 * 4 + k1] * pow(x_scaled, k1) * pow(y_scaled, k2 - 1) * float(k2) / dy;
            }
        }
    } else {
        for(int k1 = 1; k1 < 4; k1++){
            for(int k2 = 1; k2 < 4; k2++){
                result += A_vec[k2 * 4 + k1] * pow(x_scaled, k1 - 1) * pow(y_scaled, k2 - 1) * float(k1 * k2) / (dx * dy);
            }
        }
    }
    
    return result;
}

void interp::eval_all(double x, double y, double z, struct hybrid_spline_3D & spline, double & f, double df []){
    int kx, ky, kz;
    double x_eval, y_eval, z_eval, dx, dy, x_scaled, y_scaled, X_vec[16], A_vec[16];

    x_eval = in_interval(x, spline.x_vals[0], spline.x_vals[spline.length_x - 1]);  kx = find_segment(x_eval, spline.x_vals, spline.length_x, spline.accel[0]);
    y_eval = in_interval(y, spline.y_vals[0], spline.y_vals[spline.length_y - 1]);  ky = find_segment(y_eval, spline.y_vals, spline.length_y, spline.accel[1]);
    z_eval = in_interval(z, spline.z_vals[0], spline.z_vals[spline.length_z - 1]);  kz = find_segment(z_eval, spline.z_vals, spline.length_z, spline.accel[2]);
    
    dx = spline.x_vals[kx + 1] - spline.x_vals[kx], x_scaled = (x_eval - spline.x_vals[kx]) / dx;
    dy = spline.y_vals[ky + 1] - spline.y_vals[ky], y_scaled = (y_eval - spline.y_vals[ky]) / dy;
    
    // Set f, df/dx, and df/dy from bicubic of f
    X_vec[0] = eval_node_f(z_eval, spline, kx, ky, kz);                             X_vec[8] =  finite_diff_dfdy(z_eval, spline, kx, ky, kz) * dy;
    X_vec[1] = eval_node_f(z_eval, spline, kx + 1, ky, kz);                         X_vec[9] =  finite_diff_dfdy(z_eval, spline, kx + 1, ky, kz) * dy;
    X_vec[2] = eval_node_f(z_eval, spline, kx, ky + 1, kz);                         X_vec[10] = finite_diff_dfdy(z_eval, spline, kx, ky + 1, kz) * dy;
    X_vec[3] = eval_node_f(z_eval, spline, kx + 1, ky + 1, kz);                     X_vec[11] = finite_diff_dfdy(z_eval, spline, kx + 1, ky + 1, kz) * dy;
    
    X_vec[4] = finite_diff_dfdx(z_eval, spline, kx, ky, kz) * dx;                   X_vec[12] = finite_diff_ddfdxdy(z_eval, spline, kx, ky, kz) * dx * dy;
    X_vec[5] = finite_diff_dfdx(z_eval, spline, kx + 1, ky, kz) * dx;               X_vec[13] = finite_diff_ddfdxdy(z_eval, spline, kx + 1, ky, kz) * dx * dy;
    X_vec[6] = finite_diff_dfdx(z_eval, spline, kx, ky + 1, kz) * dx;               X_vec[14] = finite_diff_ddfdxdy(z_eval, spline, kx, ky + 1, kz) * dx * dy;
    X_vec[7] = finite_diff_dfdx(z_eval, spline, kx + 1, ky + 1, kz) * dx;           X_vec[15] = finite_diff_ddfdxdy(z_eval, spline, kx + 1, ky + 1, kz) * dx * dy;
    
    f = 0.0; df[0] = 0.0; df[1] = 0.0;
    for(int j = 0; j < 16; j++){
        A_vec[j] = 0.0;
        for(int k = 0; k < 16; k++){
            A_vec[j] += bicubic_conversion_matrix[j][k]*X_vec[k];
        }
    }

    for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            f += A_vec[k1 + 4 * k2] * pow(x_scaled, k1) * pow(y_scaled, k2);
        }
    }
    
    for(int k1 = 1; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            df[0] += A_vec[k1 + 4 * k2] * pow(x_scaled, k1 - 1) * pow(y_scaled, k2) * float(k1) / dx;
        }
    }

    for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 1; k2 < 4; k2++){
            df[1] += A_vec[k1 + 4 * k2] * pow(x_scaled, k1) * pow(y_scaled, k2 - 1) * float(k2) / dy;
        }
    }
	
    // Set dfdz from bicubic of dfdz
    X_vec[0] = eval_node_dfdz(z_eval, spline, kx, ky, kz);                          X_vec[8] =  finite_diff_ddfdydz(z_eval, spline, kx, ky, kz) * dy;
    X_vec[1] = eval_node_dfdz(z_eval, spline, kx + 1, ky, kz);                      X_vec[9] =  finite_diff_ddfdydz(z_eval, spline, kx + 1, ky, kz) * dy;
    X_vec[2] = eval_node_dfdz(z_eval, spline, kx, ky + 1, kz);                      X_vec[10] = finite_diff_ddfdydz(z_eval, spline, kx, ky + 1, kz) * dy;
    X_vec[3] = eval_node_dfdz(z_eval, spline, kx + 1, ky + 1, kz);                  X_vec[11] = finite_diff_ddfdydz(z_eval, spline, kx + 1, ky + 1, kz) * dy;
    
    X_vec[4] = finite_diff_ddfdxdz(z_eval, spline, kx, ky, kz) * dx;                X_vec[12] = finite_diff_dddfdxdydz(z_eval, spline, kx, ky, kz) * dx * dy;
    X_vec[5] = finite_diff_ddfdxdz(z_eval, spline, kx + 1, ky ,kz) * dx;            X_vec[13] = finite_diff_dddfdxdydz(z_eval, spline, kx + 1, ky, kz) * dx * dy;
    X_vec[6] = finite_diff_ddfdxdz(z_eval, spline, kx, ky + 1, kz) * dx;            X_vec[14] = finite_diff_dddfdxdydz(z_eval, spline, kx, ky + 1, kz) * dx * dy;
    X_vec[7] = finite_diff_ddfdxdz(z_eval, spline, kx + 1, ky + 1, kz) * dx;        X_vec[15] = finite_diff_dddfdxdydz(z_eval, spline, kx + 1, ky + 1, kz) * dx * dy;

    df[2] = 0.0;
    for(int j = 0; j < 16; j++){
        A_vec[j] = 0;
        for(int k = 0; k < 16; k++){
            A_vec[j] += bicubic_conversion_matrix[j][k] * X_vec[k];
        }
    }
    
    for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            df[2] += A_vec[k1 + 4 * k2] * pow(x_scaled, k1) * pow(y_scaled, k2);
        }
    }
}

void interp::eval_all(double x, double y, double z, struct hybrid_spline_3D & spline, double & f, double df [], double ddf []){
    int kx, ky, kz;
    double x_eval, y_eval, z_eval, dx, dy, x_scaled, y_scaled, X_vec[16], A_vec[16];
    
    x_eval = in_interval(x, spline.x_vals[0], spline.x_vals[spline.length_x - 1]);  kx = find_segment(x_eval, spline.x_vals, spline.length_x, spline.accel[0]);
    y_eval = in_interval(y, spline.y_vals[0], spline.y_vals[spline.length_y - 1]);  ky = find_segment(y_eval, spline.y_vals, spline.length_y, spline.accel[1]);
    z_eval = in_interval(z, spline.z_vals[0], spline.z_vals[spline.length_z - 1]);  kz = find_segment(z_eval, spline.z_vals, spline.length_z, spline.accel[2]);
    
    dx = spline.x_vals[kx + 1] - spline.x_vals[kx], x_scaled = (x_eval - spline.x_vals[kx]) / dx;
    dy = spline.y_vals[ky + 1] - spline.y_vals[ky], y_scaled = (y_eval - spline.y_vals[ky]) / dy;
    
    // Set f, df/dx, df/dy, and ddf/dxdy from bicubic of f
    X_vec[0] = eval_node_f(z_eval, spline, kx, ky, kz);                             X_vec[8] =  finite_diff_dfdy(z_eval, spline, kx, ky, kz) * dy;
    X_vec[1] = eval_node_f(z_eval, spline, kx + 1, ky, kz);                         X_vec[9] =  finite_diff_dfdy(z_eval, spline, kx + 1, ky, kz) * dy;
    X_vec[2] = eval_node_f(z_eval, spline, kx, ky + 1, kz);                         X_vec[10] = finite_diff_dfdy(z_eval, spline, kx, ky + 1, kz) * dy;
    X_vec[3] = eval_node_f(z_eval, spline, kx + 1, ky + 1, kz);                     X_vec[11] = finite_diff_dfdy(z_eval, spline, kx + 1, ky + 1, kz) * dy;
    
    X_vec[4] = finite_diff_dfdx(z_eval, spline, kx, ky, kz) * dx;                   X_vec[12] = finite_diff_ddfdxdy(z_eval, spline, kx, ky, kz) * dx * dy;
    X_vec[5] = finite_diff_dfdx(z_eval, spline, kx + 1, ky, kz) * dx;               X_vec[13] = finite_diff_ddfdxdy(z_eval, spline, kx + 1, ky, kz) * dx * dy;
    X_vec[6] = finite_diff_dfdx(z_eval, spline, kx, ky + 1, kz) * dx;               X_vec[14] = finite_diff_ddfdxdy(z_eval, spline, kx, ky + 1, kz) * dx * dy;
    X_vec[7] = finite_diff_dfdx(z_eval, spline, kx + 1, ky + 1, kz) * dx;           X_vec[15] = finite_diff_ddfdxdy(z_eval, spline, kx + 1, ky + 1, kz) * dx * dy;
    
    for(int j = 0; j < 16; j++){
        A_vec[j] = 0.0;
        for(int k = 0; k < 16; k++){
            A_vec[j] += bicubic_conversion_matrix[j][k]*X_vec[k];
        }
    }
    
    f = 0.0;
    for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            f += A_vec[k1 + 4 * k2] * pow(x_scaled, k1) * pow(y_scaled, k2);
        }
    }

    df[0] = 0.0;
    for(int k1 = 1; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            df[0] += A_vec[k1 + 4 * k2] * pow(x_scaled, k1 - 1) * pow(y_scaled, k2) * float(k1) / dx;
        }
    }
    
    df[1] = 0.0;
    for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 1; k2 < 4; k2++){
            df[1] += A_vec[k1 + 4 * k2] * pow(x_scaled, k1) * pow(y_scaled, k2 - 1) * float(k2) / dy;
        }
    }

    ddf[3] = 0.0;
    for(int k1 = 1; k1 < 4; k1++){
        for(int k2 = 1; k2 < 4; k2++){
            ddf[3] += A_vec[k1 + 4 * k2] * pow(x_scaled, k1 - 1) * pow(y_scaled, k2 - 1) * float(k1 * k2) / (dx * dy);
        }
    }
    
    // Set dfdz, ddf/dxdz, and ddf/dydz from bicubic of dfdz
    X_vec[0] = eval_node_dfdz(z_eval, spline, kx, ky, kz);                          X_vec[8] =  finite_diff_ddfdydz(z_eval, spline, kx, ky, kz) * dy;
    X_vec[1] = eval_node_dfdz(z_eval, spline, kx + 1, ky, kz);                      X_vec[9] =  finite_diff_ddfdydz(z_eval, spline, kx + 1, ky, kz) * dy;
    X_vec[2] = eval_node_dfdz(z_eval, spline, kx, ky + 1, kz);                      X_vec[10] = finite_diff_ddfdydz(z_eval, spline, kx, ky + 1, kz) * dy;
    X_vec[3] = eval_node_dfdz(z_eval, spline, kx + 1, ky + 1, kz);                  X_vec[11] = finite_diff_ddfdydz(z_eval, spline, kx + 1, ky + 1, kz) * dy;
    
    X_vec[4] = finite_diff_ddfdxdz(z_eval, spline, kx, ky, kz) * dx;                X_vec[12] = finite_diff_dddfdxdydz(z_eval, spline, kx, ky, kz) * dx * dy;
    X_vec[5] = finite_diff_ddfdxdz(z_eval, spline, kx + 1, ky ,kz) * dx;            X_vec[13] = finite_diff_dddfdxdydz(z_eval, spline, kx + 1, ky, kz) * dx * dy;
    X_vec[6] = finite_diff_ddfdxdz(z_eval, spline, kx, ky + 1, kz) * dx;            X_vec[14] = finite_diff_dddfdxdydz(z_eval, spline, kx, ky + 1, kz) * dx * dy;
    X_vec[7] = finite_diff_ddfdxdz(z_eval, spline, kx + 1, ky + 1, kz) * dx;        X_vec[15] = finite_diff_dddfdxdydz(z_eval, spline, kx + 1, ky + 1, kz) * dx * dy;
    
    for(int j = 0; j < 16; j++){
        A_vec[j] = 0;  
        for(int k = 0; k < 16; k++){
            A_vec[j] += bicubic_conversion_matrix[j][k] * X_vec[k];
        }
    }

    df[2] = 0.0;
    for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            df[2] += A_vec[k1 + 4 * k2] * pow(x_scaled, k1) * pow(y_scaled, k2);
        }
    }

    ddf[4] = 0.0;
    for(int k1 = 1; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            ddf[4] += A_vec[k1 + 4 * k2] * pow(x_scaled, k1 - 1) * pow(y_scaled, k2) * float(k1) / dx;
        }
    }

    ddf[5] = 0.0;
    for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 1; k2 < 4; k2++){
            ddf[5] += A_vec[k1 + 4 * k2] * pow(x_scaled, k1) * pow(y_scaled, k2 - 1) * float(k2) / dy;
        }
    }
    
    // Set ddf/dxdx from d/dx of bicubic of df/dx
    X_vec[0] = eval_node_dfdx(z_eval, spline, kx, ky, kz);                          X_vec[8] =  finite_diff_ddfdxdy(z_eval, spline, kx, ky, kz) * dy;
    X_vec[1] = eval_node_dfdx(z_eval, spline, kx + 1, ky, kz);                      X_vec[9] =  finite_diff_ddfdxdy(z_eval, spline, kx + 1, ky, kz) * dy;
    X_vec[2] = eval_node_dfdx(z_eval, spline, kx, ky + 1, kz);                       X_vec[10] = finite_diff_ddfdxdy(z_eval, spline, kx, ky + 1, kz) * dy;
    X_vec[3] = eval_node_dfdx(z_eval, spline, kx + 1, ky + 1, kz);                   X_vec[11] = finite_diff_ddfdxdy(z_eval, spline, kx + 1, ky + 1, kz) * dy;
    
    X_vec[4] = finite_diff_ddfdxdx(z_eval, spline, kx, ky, kz) * dx;                X_vec[12] = finite_diff_dddfdxdxdy(z_eval, spline, kx, ky, kz) * dx * dy;
    X_vec[5] = finite_diff_ddfdxdx(z_eval, spline, kx + 1, ky, kz) * dx;            X_vec[13] = finite_diff_dddfdxdxdy(z_eval, spline, kx + 1, ky, kz) * dx * dy;
    X_vec[6] = finite_diff_ddfdxdx(z_eval, spline, kx, ky + 1, kz) * dx;            X_vec[14] = finite_diff_dddfdxdxdy(z_eval, spline, kx, ky + 1, kz) * dx * dy;
    X_vec[7] = finite_diff_ddfdxdx(z_eval, spline, kx + 1, ky + 1, kz) * dx;        X_vec[15] = finite_diff_dddfdxdxdy(z_eval, spline, kx + 1, ky + 1, kz) * dx * dy;
    
    for(int j = 0; j < 16; j++){
        A_vec[j] = 0;
        for(int k = 0; k < 16; k++){
            A_vec[j] += bicubic_conversion_matrix[j][k]*X_vec[k];
        }
    }

    ddf[0] = 0.0;
    for(int k1 = 1; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            ddf[0] += A_vec[k1 + 4 * k2] * pow(x_scaled, k1 - 1) * pow(y_scaled, k2) * float(k1) / dx;
        }
    }
    
    // Set ddf/dydy from d/dy of bicubic of df/dy
    X_vec[0] = eval_node_dfdy(z_eval, spline, kx, ky, kz);                          X_vec[8] =  finite_diff_ddfdydy(z_eval, spline, kx, ky, kz) * dy;
    X_vec[1] = eval_node_dfdy(z_eval, spline, kx + 1, ky, kz);                      X_vec[9] =  finite_diff_ddfdydy(z_eval, spline, kx + 1, ky, kz) * dy;
    X_vec[2] = eval_node_dfdy(z_eval, spline, kx, ky + 1, kz);                      X_vec[10] = finite_diff_ddfdydy(z_eval, spline, kx, ky + 1, kz) * dy;
    X_vec[3] = eval_node_dfdy(z_eval, spline, kx + 1, ky + 1, kz);                  X_vec[11] = finite_diff_ddfdydy(z_eval, spline, kx + 1, ky + 1, kz) * dy;
    
    X_vec[4] = finite_diff_ddfdxdy(z_eval, spline, kx, ky, kz) * dx;                X_vec[12] = finite_diff_dddfdxdydy(z_eval, spline, kx, ky, kz) * dx * dy;
    X_vec[5] = finite_diff_ddfdxdy(z_eval, spline, kx + 1, ky, kz) * dx;            X_vec[13] = finite_diff_dddfdxdydy(z_eval, spline, kx + 1, ky, kz) * dx * dy;
    X_vec[6] = finite_diff_ddfdxdy(z_eval, spline, kx, ky + 1, kz) * dx;            X_vec[14] = finite_diff_dddfdxdydy(z_eval, spline, kx, ky + 1, kz) * dx * dy;
    X_vec[7] = finite_diff_ddfdxdy(z_eval, spline, kx + 1, ky + 1, kz) * dx;        X_vec[15] = finite_diff_dddfdxdydy(z_eval, spline, kx + 1, ky + 1, kz) * dx * dy;
    
    for(int j = 0; j < 16; j++){
        A_vec[j] = 0;
        for(int k = 0; k < 16; k++){
            A_vec[j] += bicubic_conversion_matrix[j][k]*X_vec[k];
        }
    }

    ddf[1] = 0.0;
    for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 1; k2 < 4; k2++){
            ddf[1] += A_vec[k1 + 4 * k2] * pow(x_scaled, k1) * pow(y_scaled, k2 - 1) * float(k2) / dy;
        }
    }
    
    // Set ddf/dzdz from bicubic of ddf/dzdz
    X_vec[0] = eval_node_ddfdzdz(z, spline, kx, ky, kz);                            X_vec[8] =  finite_diff_dddfdydzdz(z, spline, kx, ky, kz) * dy;
    X_vec[1] = eval_node_ddfdzdz(z, spline, kx + 1, ky, kz);                        X_vec[9] =  finite_diff_dddfdydzdz(z, spline, kx + 1, ky, kz) * dy;
    X_vec[2] = eval_node_ddfdzdz(z, spline, kx, ky + 1, kz);                        X_vec[10] = finite_diff_dddfdydzdz(z, spline, kx, ky + 1, kz) * dy;
    X_vec[3] = eval_node_ddfdzdz(z, spline, kx + 1, ky + 1, kz);                    X_vec[11] = finite_diff_dddfdydzdz(z, spline, kx + 1, ky + 1, kz) * dy;
    
    X_vec[4] = finite_diff_dddfdxdzdz(z, spline, kx, ky, kz) * dx;                  X_vec[12] = finite_diff_ddddfdxdydzdz(z, spline, kx, ky, kz) * dx * dy;
    X_vec[5] = finite_diff_dddfdxdzdz(z, spline, kx + 1, ky, kz) * dx;              X_vec[13] = finite_diff_ddddfdxdydzdz(z, spline, kx + 1, ky, kz) * dx * dy;
    X_vec[6] = finite_diff_dddfdxdzdz(z, spline, kx, ky + 1, kz) * dx;              X_vec[14] = finite_diff_ddddfdxdydzdz(z, spline, kx, ky + 1, kz) * dx * dy;
    X_vec[7] = finite_diff_dddfdxdzdz(z, spline, kx + 1, ky + 1, kz) * dx;          X_vec[15] = finite_diff_ddddfdxdydzdz(z, spline, kx + 1, ky + 1, kz) * dx * dy;

    for(int j = 0; j < 16; j++){
        A_vec[j] = 0;
        for(int k = 0; k < 16; k++){
            A_vec[j] += bicubic_conversion_matrix[j][k]*X_vec[k];
        }
    }

    ddf[2] = 0;
    for(int k1 = 0; k1 < 4; k1++){
        for(int k2 = 0; k2 < 4; k2++){
            ddf[2] += A_vec[k1 + 4 * k2] * pow(x_scaled, k1) * pow(y_scaled, k2);
        }
    }
}




#endif /* _INTERPOLATION_CPP_ */
