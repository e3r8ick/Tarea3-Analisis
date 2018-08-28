/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 04.08.2018
 */

#include <cmath>
#include <math.h>
#include <limits>
#include <functional>

#include "Exception.hpp"

#ifndef ANPI_ROOT_RIDDER_HPP
#define ANPI_ROOT_RIDDER_HPP

namespace anpi {

  /**
   * @param a 
   * @param b
   *
   * @return a number with the same magnitude as "a" and the same sign as "b".
   */
  template<typename T>
  T Sign(T a, T b){
    T temp = 0;
    if(signbit(b)){
      temp = a*(-1);
    }
    else{
      temp = std::abs(a);
    }
    return temp;
  }
  
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
  T rootRidder(const std::function<T(T)>& funct,T xi,T xii,const T eps) {
    const int MAXIT=60;
    T fl=funct(xi);
    T fh=funct(xii);
    if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
      T xl=xi;
      T xh=xii;
      T ans=-9.99e99;
      for (int j = 0; j < MAXIT; j++) {
        T xm=0.5*(xl+xh);
        T fm=funct(xm);
        T s=std::sqrt(fm*fm-fl*fh);
        if (s == 0.0) {
          return ans;
        }
        T xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s); 
        if (std::abs(xnew-ans) <= eps){ 
          return ans;
        }
        ans=xnew;
        T fnew=funct(ans);
        if (fnew == 0.0) {
          return ans;
        }
        if (anpi::Sign(fm,fnew) != fm) {
          xl=xm;
          fl=fm;
          xh=ans;
          fh=fnew;
        } else if (anpi::Sign(fl,fnew) != fl) {
          xh=ans;
          fh=fnew;
        } else if (anpi::Sign(fh,fnew) != fh) {
          xl=ans;
          fl=fnew;
        } else {
          throw("never get here.");
        }
        if (std::abs(xh-xl) <= eps) {
          std::cout << "error8: "<< std::endl;
          return ans;
        }
      }
      throw("rootRidder exceed maximum iterations");
    }
    else {
      std::cout << "error9: "<< std::endl;
      if (fl == 0.0) {
        std::cout << "error10: "<< std::endl;
        return xi;
      }
      if (fh == 0.0) {
        std::cout << "error11: "<< std::endl;
        return xii;
      }
      throw("root must be bracketed in rootRidder.");
    }
    return std::numeric_limits<T>::quiet_NaN();
  }

}
  
#endif

