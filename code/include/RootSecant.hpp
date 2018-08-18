/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 */

#include <cmath>
#include <limits>
#include <functional>

#include "Exception.hpp"

#ifndef ANPI_ROOT_SECANT_HPP
#define ANPI_ROOT_SECANT_HPP

namespace anpi {
  
  /**
   * Find a root of the function funct looking for it starting at xi
   * by means of the secant method.
   *
   * @param funct a functor of the form "T funct(T x)"
   * @param xi initial position
   * @param xii second initial position 
   *
   * @return root found, or NaN if no root could be found
   */
  template<typename T>
  T rootSecant(const std::function<T(T)>& funct,T xi,T xii,const T eps) {
    T fl = funct(xii);
    T f = funct(xi);
    int maxi = std::numeric_limits<T>::digits;
    T ea(T(0));
    for (int j = 0; j < maxi; ++j){
      T dx = (xii - xi)*f/(f-fl);
      xii = xi;
      fl = f;
      xi += dx;
      f = funct(xi);
      if (std::abs(xi) > eps){
        ea = std::abs((xi-xii)/xi)*T(100); 
      }
      if (ea < eps){
        return xi;
      }
    }
    // Return NaN if no root was found
    return std::numeric_limits<T>::quiet_NaN();
  }

}
  
#endif

