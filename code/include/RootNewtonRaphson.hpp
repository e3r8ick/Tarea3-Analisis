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

#ifndef ANPI_NEWTON_RAPHSON_HPP
#define ANPI_NEWTON_RAPHSON_HPP

namespace anpi {
  
  /**
   * Find the roots of the function funct looking by means of the
   * Newton-Raphson method
   *
   * @param funct a functor of the form "T funct(T x)"
   * @param xi initial root guess
   * 
   * @return root found, or NaN if none could be found.
   *
   * @throws anpi::Exception if inteval is reversed or both extremes
   *         have same sign.
   */

  template<typename T>
  T rootNewtonRaphson(const std::function<T(T)>& funct,T xi,const T eps) {
    T f = funct(xi);
    T h = T(0.1);
    T df = (funct(xi+h) - funct(xi-h))/(2*h);
    T maxi = std::numeric_limits<T>::digits;
    maxi = pow(maxi,3);
    T ea(T(0));
    for (T j = maxi; j > 0; --j){
      T xiold = xi;
      if (std::abs(df) < eps){
        throw anpi::Exception("Division sobre 0");
      }
      T divi = f/df;
      xi = xi - divi;
      if (std::abs(xi) > eps){
        ea = std::abs((xi-xiold)/xi)*T(100); 
      }
      f = funct(xi);
      df = (funct(xi+h) - funct(xi-h))/(2*h);
      if (ea < std::sqrt(eps) && std::abs(f) < eps){ 
        return xi;
      }
    }
    // Return NaN if no root was found
    std::cout << "Newton Raphson: no root was found for " << eps << std::endl;
    return std::numeric_limits<T>::quiet_NaN();
  }

}
  
#endif
