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

#ifndef ANPI_ROOT_BRENT_HPP
#define ANPI_ROOT_BRENT_HPP

namespace anpi {

  /**
   * Find the roots of the function funct looking for it in the
   * interval [xl,xu], using the Brent's method.
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
  T rootBrent(const std::function<T(T)>& funct,T xl,T xu,const T eps) {

    T fl = funct(xl);
    T fu = funct(xu);
    if (xl > xu){
      throw anpi::Exception("Interval reversed");
    }
    if (std::signbit(fl) == std::signbit(fu)){
      throw anpi::Exception("Signos iguales");
    }
    T root;
    int maxi = std::numeric_limits<T>::digits;
    T x1 = xl;
    T x2 = (xl+xu)/2;
    T x3 = xu;
    T y1 = funct(x1);
    T y2 = funct(x2);
    T y3 = funct(x3);
    T xN;
    T ea (T(0));
    maxi = maxi*maxi;
    for (int j = maxi; j > 0; --j){
      xN = x1*y2*y3/((y1-y2)*(y1-y3)) + x2*y1*y3/((y2-y1)*(y2-y3)) + x3*y1*y2/((y3-y1)*(y3-y2));
      if (std::abs(xN) > eps){
        ea = std::abs((xN-x3)/xN)*T(100);
      }
      y1 = y2;
      y2 = y3;
      y3 = funct(xN);
      x1 = x2;
      x2 = x3;
      x3 = xN;

      if (ea < eps && (std::abs(y3) - eps) < 0){
        return x3;
      }
    }    
    root = anpi::rootSecant(funct, xl, xu, eps);
    if (!std::isnan(root)){
      return root;
    }
    root = anpi::rootBisection(funct, xl, xu, eps);
    if (!std::isnan(root)){
      std::cout << "Brent: Converge por Bi" << std::endl;
      return root;
    }
    // Return NaN if no root was found
    return std::numeric_limits<T>::quiet_NaN();
  }
}
  
#endif

