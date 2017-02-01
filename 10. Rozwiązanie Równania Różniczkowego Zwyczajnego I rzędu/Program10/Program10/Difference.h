#ifndef __DIFFERENCE
#define __DIFFERENCE

#include<iostream>

class DifferenceException : public std::runtime_error {
	public:
		explicit DifferenceException(const std::string& str) : std::runtime_error(str) {}
};

/** Difference - class helps to calculate difference - approximation of derivatives.
* @author Christopher Abram
*/
template<typename T>
class Difference {
		
		// Aliases:
		using value_type = T;
		using self = Difference<T>;
		using self_reference = self&;

		// Increase:
		value_type __h;

		// Real function of one real variable:
		typedef value_type(*function)(value_type);
		function __f;

	public:

		using value_type = T;

		// Specifies special argument, to decide which difference calculate:
		enum Points {
			TWO,
			THREE
		};

		/** Calculates forward difference.
		* @param	value_type x		- the point to calculate derivative of function f in.
		* @param	Points special		- special flag to determine which method to choose, two, three or more points...
		* @return	value_type			- the difference in given point for given function.
		* @throws	DifferenceException - when increase is default (zero).
		*/
		value_type forward(value_type x, Points special = Points::TWO){
			value_type diff{};
			function f = getFunction();
			value_type h = getIncrease();

			if (h == value_type{})
				throw *(new DifferenceException("Error: Unable to calculate forward difference, because increase value is default (zero)"));

			/* Two-point difference for forward method */
			if (special == Points::TWO)
				diff = (f(x + h) - f(x)) / h;
			/* Three-point difference for forward method */
			else if (special == Points::THREE)
				diff = ((value_type)(-1.5) * f(x) + (value_type)2 * f(x + h) + (value_type)(-0.5) * f(x + h + h)) / h;

			return diff;
		}// end forward


		 /** Calculates central difference.
		 * @param	value_type x		- the point to calculate derivative of function f in.
		 * @param	Points special		- special flag to determine which method to choose, two, three or more points...
		 * @return	value_type			- the difference in given point for given function.
		 * @throws	DifferenceException - when increase is default (zero).
		 */
		value_type central(value_type x, Points special = Points::TWO) {
			value_type diff{};
			function f = getFunction();
			value_type h = getIncrease();

			if (h == value_type{})
				throw *(new DifferenceException("Error: Unable to calculate central difference, because increase value is default (zero)"));

			/* Two-point difference for central method */
			if (special == Points::TWO)
				diff = (f(x + h) - f(x - h)) / (h + h);
			/* Three-point difference for central method */
			else if (special == Points::THREE) {
				// TODO: find approximation ...
			}
			return diff;
		}// end central


		 /** Calculates backward difference.
		 * @param	value_type x		- the point to calculate derivative of function f in.
		 * @param	Points special		- special flag to determine which method to choose, two, three or more points...
		 * @return	value_type			- the difference in given point for given function.
		 * @throws	DifferenceException - when increase is default (zero).
		 */
		value_type backward(value_type x, Points special = Points::TWO) {
			value_type diff{};
			function f = getFunction();
			value_type h = getIncrease();

			if (h == value_type{})
				throw *(new DifferenceException("Error: Unable to calculate backward difference, because increase value is default (zero)"));

			/* Two-point difference for backward method */
			if (special == Points::TWO)
				diff = (f(x) - f(x - h)) / h;
			/* Three-point difference for backward method */
			else if (special == Points::THREE)
				diff = ((value_type)(0.5) * f(x - h - h) - (value_type)2 * f(x - h) + (value_type)(1.5) * f(x)) / h;

			return diff;
		}// end backward






		self_reference setFunction(function f) {
			__f = f;
			return *this;
		}// end setFunction

		function getFunction(void) const {
			return __f;
		}// end getFunction

		self_reference setIncrease(value_type h) {
			if (h == value_type{})
				throw *(new DifferenceException("Error: Increase value can not be default (zero)."));
			__h = h;
			return *this;
		}// end setIncrease

		value_type getIncrease(void) const {
			return __h;
		}// end getIncrease
};

#endif /* __DIFFERENCE */