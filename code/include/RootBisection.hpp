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
#include <iostream>

#include "Exception.hpp"

#ifndef ANPI_ROOT_BISECTION_HPP
#define ANPI_ROOT_BISECTION_HPP

namespace anpi {
  
  /**
   * Find the roots of the function funct looking for it in the
   * interval [xl,xu], using the bisection method.
   *
   * @param funct a std::function of the form "T funct(T x)"
   * @param xl lower interval limit
   * @param xu upper interval limit
   *
   * @return root found, or NaN if none could be found.
   *
   * @throws anpi::Exception if inteval is reversed or both extremes
   *         have same sign.
   */
  template<typename T>
  T rootBisection(const std::function<T(T)>& funct,T xl,T xu,const T eps) {

    T fl = funct(xl);
    T fu = funct(xu);
    if (xl > xu){
      throw anpi::Exception("Interval reversed");
    }
    if (std::signbit(fl) == std::signbit(fu)){
      throw anpi::Exception("Signos iguales");
    }
    T xr = xl;
    T fr;
    T ea = T();
    T xold;
    int maxi = std::numeric_limits<T>::digits;
    maxi *= maxi;
    for (int i = maxi; i > 0 ; --i){
      xold = xr;
      xr = (xl+xu)/2;
      fr = funct(xr);
      T cond = fl*fr;
      if (std::abs(xr) > eps){ //evita div por 0
        ea = std::abs((xr-xold)/xr)*T(100); //error aprox
      }
      if (cond < T(0)){
        xu = xr;
      }else if (cond > T(0)){
        xl = xr;
        fl = fr;
      }else{
        ea = T(0);
        xr = (std::abs(fl) < eps) ? xl : xr; // fl == 0
      }
      if (ea < eps){
        return xr;  // retorna si se obtiene la precision
      }
    }
    // Return NaN if no root was found
    return std::numeric_limits<T>::quiet_NaN();
  }

}
  
#endif

