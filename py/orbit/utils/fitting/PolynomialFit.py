# -------------------------------------------------------------------------
# This is a polynom fitting class. It uses the standart
# LSM nmatrix approach.
# As input you can use Function or SplineCH instances
# order of the polynomial a0+...+aN*x^N is N
# --------------------------------------------------------------------------
import sys
import math

from orbit.core.orbit_utils import Function, SplineCH
from orbit.core.orbit_utils import PhaseVector, Matrix
from orbit.core.orbit_utils import Polynomial


class PolynomialFit:
    def __init__(self, order):
        """The calss to fit Function or SplineCH with a polynomial."""
        self.order = order
        self.polynomial = Polynomial()
        # self.coef_err_arr is a final array with coef_arr and err_arr for polinomial coefficients
        self.coef_err_arr = []
        # self.x_y_err_arr is initial data with (x,y,y_err) points
        self.x_y_err_arr = []

    def getPolynomial(self):
        """It returns an unbounded polynomial."""
        polynomial = Polynomial()
        self.polynomial.copyTo(polynomial)
        return polynomial

    def getCoefficientsAndErr(self):
        """It returns the list of coefficients and their errors"""
        return self.coef_err_arr

    def getCoefficients(self):
        """Returns the list of coefficients of the polynomial"""
        coef_arr = []
        for i in range(len(self.coef_err_arr)):
            [coef, err] = self.coef_err_arr[i]
            coef_arr.append(coef)
        return coef_arr

    def fitFunction(self, function):
        """Fit the Function instance"""
        self.x_y_err_arr = []
        for i in range(function.getSize()):
            x = function.x(i)
            y = function.y(i)
            err = function.err(i)
            self.x_y_err_arr.append([x, y, err])
        self._makePolynomial()

    def fitSpline(self, spline):
        """Fit the SplineCH instance"""
        self.x_y_err_arr = []
        for i in range(spline.getSize()):
            x = spline.x(i)
            y = spline.y(i)
            err = 0.0
            self.x_y_err_arr.append([x, y, err])
        self._makePolynomial()

    def _makePolynomial(self):
        nPoints = len(self.x_y_err_arr)
        if(nPoints < (self.order + 1)):
            self.order = nPoints - 1
        #---- now make A and A^T matrices and y-vector
        aMtr = Matrix(nPoints, self.order + 1)
        yVct = PhaseVector(nPoints)
        for i in range(nPoints):
            y = self.x_y_err_arr[i][1]
            yVct.set(i,y)
            for j in range(self.order + 1):
                x = self.x_y_err_arr[i][0]
                aMtr.set(i, j, math.pow(x, j))
        aMtrT = aMtr.transpose()
        # now the resuting coefficients and errors
        ATA_inv = (aMtrT.mult(aMtr)).invert()
        coeffVct = (ATA_inv.mult(aMtrT)).mult(yVct)
        # polinimial coefficients have been found -> coeffVct
        coef_arr = [0.0] * (self.order + 1)   
        self.polynomial.order(self.order)
        for i in range(self.order + 1):
            coef_arr[i] = coeffVct.get(i)
            self.polynomial.coefficient(i,coef_arr[i])
        #---- here we estimate errors by deviation y_fit from y_init only
        #---- It means that we do not pay attention to different weights
        #---- of the initial Y data errors
        coeff_err_arr = [0.0] * (self.order + 1)
        for i in range(self.order + 1):
           coeff_err_arr[i] = math.sqrt(ATA_inv.get(i,i))
        total_sigma = 0.0
        for k in range(nPoints):
            x = self.x_y_err_arr[k][0]
            y = self.x_y_err_arr[k][1]
            total_sigma += (self.polynomial.value(x) - y) ** 2
        degrees_of_freedom = int(abs(nPoints - (self.order + 1)))
        if(degrees_of_freedom == 0): degrees_of_freedom = 1
        total_sigma = math.sqrt(total_sigma / degrees_of_freedom)
        for i in range(self.order + 1):
            coeff_err_arr[i] *= total_sigma
        # set the resulting coefficients and errors array
        self.coef_err_arr = [coef_arr, coeff_err_arr]
