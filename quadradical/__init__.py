__all__ = ["QuadraticRadical"]

import math
import sys
from fractions import Fraction
from numbers import Integral, Rational, Real, Complex
from typing import Dict, List, NoReturn, Tuple, Union, Optional, overload

if sys.version_info >= (3, 8):
    from typing import Literal

_Complex = Union[Complex, complex]
_Real = Union[Real, float, int]
_Rational = Union[Rational, int]
_Integral = Union[Integral, int]


def _sqrt(x: int, times: int = 5, initial: Optional[Fraction] = None) \
        -> Fraction:
    res: Fraction = Fraction(x >> (x.bit_length()//2)) if initial is None \
        else initial
    for _ in range(times):
        res -= (res ** 2 - x) / (res * 2)
        a = int(res)
        if a ** 2 == x:
            return Fraction(a)
    return res


def _simplify(root: int) -> Tuple[int, int]:
    resl: int = 1
    resr: int = root
    for i in range(2, int(_sqrt(root)) + 1):
        while not resr % i ** 2:
            resl *= i
            resr //= i ** 2
        if resr == 1:
            break
    return resl, resr


class QuadraticRadical(Real):
    """Quadratic Radical class"""

    __slot__ = ["quantity"]

    quantity: Dict[int, Fraction]

    def __init__(
            self, *args: Union[_Rational, Tuple[_Rational, _Integral]]
    ) -> None:
        """"""

        self.quantity = {}
        for i in args:
            if isinstance(i, tuple):
                simplify_result: Tuple[int, int] = _simplify(int(i[1]))
                if simplify_result[1] in self.quantity:
                    self.quantity[simplify_result[1]] += \
                        Fraction(i[0] * simplify_result[0])
                else:
                    self.quantity[simplify_result[1]] = \
                        Fraction(i[0] * simplify_result[0])
            elif 1 in self.quantity:
                self.quantity[1] += Fraction(i)
            else:
                self.quantity[1] = Fraction(i)
        for k in list(self.quantity):
            if not self.quantity[k]:
                del self.quantity[k]

    def __repr__(self) -> str:
        rl: List[str] = []
        for k, v in self.quantity.items():
            rl.append(str((v, k)))
        return type(self).__name__ + "(" + ", ".join(rl) + ")"

    def __float__(self) -> float:
        res: float = 0.
        for k, v in self.quantity.items():
            res += float(v) * k ** 0.5
        return res

    def __int__(self) -> int:
        a: Optional[Fraction] = self.rational()
        if a is not None:
            return int(a)
        sqrt_approx_ge: Dict[int, Fraction] = {}
        sqrt_approx_lt: Dict[int, Fraction] = {}
        for i in self.quantity:
            sqrt_approx_ge[i] = _sqrt(i)
            if sqrt_approx_ge[i].as_integer_ratio()[1] == 1:
                sqrt_approx_lt[i] = sqrt_approx_ge[i]
            else:
                a = _sqrt(i, 1, sqrt_approx_ge[i])
                sqrt_approx_lt[i] = a * 2 - sqrt_approx_ge[i]
                sqrt_approx_ge[i] = a
        while True:
            sum_approx_lt = sum((v*sqrt_approx_ge[k] if v < 0 else
                                 v*sqrt_approx_lt[k] for k, v in
                                 self.quantity.items()), Fraction(0))
            sum_approx_ge = sum((v*sqrt_approx_lt[k] if v < 0 else
                                 v*sqrt_approx_ge[k] for k, v in
                                 self.quantity.items()), Fraction(0))
            b = int(sum_approx_lt)
            if b == int(sum_approx_ge):
                return b
            for i in self.quantity:
                if sqrt_approx_ge[i].as_integer_ratio()[1] != 1:
                    a = _sqrt(i, 1, sqrt_approx_ge[i])
                    sqrt_approx_lt[i] = a if a.as_integer_ratio()[1] == 1 \
                        else a * 2 - sqrt_approx_ge[i]
                    sqrt_approx_ge[i] = a

    def __floor__(self) -> int:
        a: Optional[Fraction] = self.rational()
        if a is not None:
            return math.floor(a)
        sqrt_approx_ge: Dict[int, Fraction] = {}
        sqrt_approx_lt: Dict[int, Fraction] = {}
        for i in self.quantity:
            sqrt_approx_ge[i] = _sqrt(i)
            if sqrt_approx_ge[i].as_integer_ratio()[1] == 1:
                sqrt_approx_lt[i] = sqrt_approx_ge[i]
            else:
                a = _sqrt(i, 1, sqrt_approx_ge[i])
                sqrt_approx_lt[i] = a * 2 - sqrt_approx_ge[i]
                sqrt_approx_ge[i] = a
        while True:
            sum_approx_lt = sum((v*sqrt_approx_ge[k] if v < 0 else
                                 v*sqrt_approx_lt[k] for k, v in
                                 self.quantity.items()), Fraction(0))
            sum_approx_ge = sum((v*sqrt_approx_lt[k] if v < 0 else
                                 v*sqrt_approx_ge[k] for k, v in
                                 self.quantity.items()), Fraction(0))
            b = math.floor(sum_approx_lt)
            if b == math.floor(sum_approx_ge):
                return b
            for i in self.quantity:
                if sqrt_approx_ge[i].as_integer_ratio()[1] != 1:
                    a = _sqrt(i, 1, sqrt_approx_ge[i])
                    sqrt_approx_lt[i] = a if a.as_integer_ratio()[1] == 1 \
                        else a * 2 - sqrt_approx_ge[i]
                    sqrt_approx_ge[i] = a

    def __ceil__(self) -> int:
        a: Optional[Fraction] = self.rational()
        if a is not None:
            return math.ceil(a)
        sqrt_approx_ge: Dict[int, Fraction] = {}
        sqrt_approx_lt: Dict[int, Fraction] = {}
        for i in self.quantity:
            sqrt_approx_ge[i] = _sqrt(i)
            if sqrt_approx_ge[i].as_integer_ratio()[1] == 1:
                sqrt_approx_lt[i] = sqrt_approx_ge[i]
            else:
                a = _sqrt(i, 1, sqrt_approx_ge[i])
                sqrt_approx_lt[i] = a * 2 - sqrt_approx_ge[i]
                sqrt_approx_ge[i] = a
        while True:
            sum_approx_lt = sum((v*sqrt_approx_ge[k] if v < 0 else
                                 v*sqrt_approx_lt[k] for k, v in
                                 self.quantity.items()), Fraction(0))
            sum_approx_ge = sum((v*sqrt_approx_lt[k] if v < 0 else
                                 v*sqrt_approx_ge[k] for k, v in
                                 self.quantity.items()), Fraction(0))
            b = math.ceil(sum_approx_lt)
            if b == math.ceil(sum_approx_ge):
                return b
            for i in self.quantity:
                if sqrt_approx_ge[i].as_integer_ratio()[1] != 1:
                    a = _sqrt(i, 1, sqrt_approx_ge[i])
                    sqrt_approx_lt[i] = a if a.as_integer_ratio()[1] == 1 \
                        else a * 2 - sqrt_approx_ge[i]
                    sqrt_approx_ge[i] = a

    def __trunc__(self) -> int:
        return int(self)

    @overload
    def __round__(self, ndigits: None = None) -> int:
        ...
    @overload
    def __round__(self, ndigits: int) -> "QuadraticRadical":
        ...
    def __round__(self, ndigits=None):
        a: Optional[Fraction] = self.rational()
        if a is not None:
            b = round(a, ndigits)
            return b if ndigits is None else QuadraticRadical(b)
        sqrt_approx_ge: Dict[int, Fraction] = {}
        sqrt_approx_lt: Dict[int, Fraction] = {}
        for i in self.quantity:
            sqrt_approx_ge[i] = _sqrt(i)
            if sqrt_approx_ge[i].as_integer_ratio()[1] == 1:
                sqrt_approx_lt[i] = sqrt_approx_ge[i]
            else:
                a = _sqrt(i, 1, sqrt_approx_ge[i])
                sqrt_approx_lt[i] = a * 2 - sqrt_approx_ge[i]
                sqrt_approx_ge[i] = a
        while True:
            sum_approx_lt = sum((v*sqrt_approx_ge[k] if v < 0 else
                                 v*sqrt_approx_lt[k] for k, v in
                                 self.quantity.items()), Fraction(0))
            sum_approx_ge = sum((v*sqrt_approx_lt[k] if v < 0 else
                                 v*sqrt_approx_ge[k] for k, v in
                                 self.quantity.items()), Fraction(0))
            b = round(sum_approx_lt, ndigits)
            if b == round(sum_approx_ge, ndigits):
                return b if ndigits is None else QuadraticRadical(b)
            for i in self.quantity:
                if sqrt_approx_ge[i].as_integer_ratio()[1] != 1:
                    a = _sqrt(i, 1, sqrt_approx_ge[i])
                    sqrt_approx_lt[i] = a if a.as_integer_ratio()[1] == 1 \
                        else a * 2 - sqrt_approx_ge[i]
                    sqrt_approx_ge[i] = a

    def __bool__(self) -> bool:
        return not not self.quantity

    def __complex__(self) -> complex:
        return float(self) + 0j

    def __hash__(self) -> int:
        return hash(self.quantity)

    def as_integer_ratio(self) -> Union[Tuple[int, int], NoReturn]:
        if not self.quantity:
            return 0, 1
        if len(self.quantity) == 1 and next(iter(self.quantity)) == 1:
            return self.quantity[1].as_integer_ratio()
        raise ValueError("An irrational cannot be written as integer ratio")

    def rational(self) -> Optional[Fraction]:
        if not self.quantity:
            return Fraction(0)
        if len(self.quantity) == 1 and next(iter(self.quantity)) == 1:
            return self.quantity[1]
        return None

    def __pos__(self) -> "QuadraticRadical":
        return self

    def __neg__(self) -> "QuadraticRadical":
        return QuadraticRadical(*((-v, k) for k, v in self.quantity.items()))

    def __abs__(self) -> "QuadraticRadical":
        return -self if self.sign() < 0 else self

    if sys.version_info >= (3, 8):
        def sign(self) -> Literal[-1, 0, 1]:
            a: Optional[Fraction] = self.rational()
            if a is not None:
                return (1 if a > 0 else -1) if a else 0
            sqrt_approx_ge: Dict[int, Fraction] = {}
            sqrt_approx_lt: Dict[int, Fraction] = {}
            for i in self.quantity:
                sqrt_approx_ge[i] = _sqrt(i)
                if sqrt_approx_ge[i].as_integer_ratio()[1] == 1:
                    sqrt_approx_lt[i] = sqrt_approx_ge[i]
                else:
                    a = _sqrt(i, 1, sqrt_approx_ge[i])
                    sqrt_approx_lt[i] = a * 2 - sqrt_approx_ge[i]
                    sqrt_approx_ge[i] = a
            while True:
                sum_approx_lt = sum((v*sqrt_approx_ge[k] if v < 0 else
                                     v*sqrt_approx_lt[k] for k, v in
                                     self.quantity.items()), Fraction(0))
                sum_approx_ge = sum((v*sqrt_approx_lt[k] if v < 0 else
                                     v*sqrt_approx_ge[k] for k, v in
                                     self.quantity.items()), Fraction(0))
                if sum_approx_lt > 0 and sum_approx_ge > 0:
                    return 1
                if sum_approx_lt < 0 and sum_approx_ge < 0:
                    return -1
                for i in self.quantity:
                    if sqrt_approx_ge[i].as_integer_ratio()[1] != 1:
                        a = _sqrt(i, 1, sqrt_approx_ge[i])
                        sqrt_approx_lt[i] = a if a.as_integer_ratio()[1] == 1 \
                            else a * 2 - sqrt_approx_ge[i]
                        sqrt_approx_ge[i] = a
    else:
        def sign(self) -> int:
            a: Optional[Fraction] = self.rational()
            if a is not None:
                return (1 if a > 0 else -1) if a else 0
            sqrt_approx_ge: Dict[int, Fraction] = {}
            sqrt_approx_lt: Dict[int, Fraction] = {}
            for i in self.quantity:
                sqrt_approx_ge[i] = _sqrt(i)
                if sqrt_approx_ge[i].as_integer_ratio()[1] == 1:
                    sqrt_approx_lt[i] = sqrt_approx_ge[i]
                else:
                    a = _sqrt(i, 1, sqrt_approx_ge[i])
                    sqrt_approx_lt[i] = a * 2 - sqrt_approx_ge[i]
                    sqrt_approx_ge[i] = a
            while True:
                sum_approx_lt = sum((-v*sqrt_approx_ge[k] if v < 0 else
                                     v*sqrt_approx_lt[k] for k, v in
                                     self.quantity.items()), Fraction(0))
                sum_approx_ge = sum((-v*sqrt_approx_lt[k] if v < 0 else
                                     v*sqrt_approx_ge[k] for k, v in
                                     self.quantity.items()), Fraction(0))
                if sum_approx_lt > 0 and sum_approx_ge > 0:
                    return 1
                if sum_approx_lt < 0 and sum_approx_ge < 0:
                    return -1
                for i in self.quantity:
                    if sqrt_approx_ge[i].as_integer_ratio()[1] != 1:
                        a = _sqrt(i, 1, sqrt_approx_ge[i])
                        sqrt_approx_lt[i] = a if a.as_integer_ratio()[1] == 1 \
                            else a * 2 - sqrt_approx_ge[i]
                        sqrt_approx_ge[i] = a

    @overload
    def __add__(self, other: Union[
        _Integral, _Rational, float, "QuadraticRadical"
    ]) -> "QuadraticRadical":
        ...
    @overload
    def __add__(self, other: Union[_Real, "QuadraticRadical"]) -> float:
        ...
    @overload
    def __add__(self, other: Complex) -> complex:
        ...
    def __add__(self, other):
        if isinstance(other, (Rational, float)):
            return self + QuadraticRadical(Fraction(*other.as_integer_ratio()))
        if not isinstance(other, QuadraticRadical):
            if isinstance(other, Real):
                return float(self) + other
            if isinstance(other, Complex):
                return complex(self) + other
            return NotImplemented
        a: Dict[int, Fraction] = self.quantity.copy()
        for k, v in other.quantity.items():
            if k in a:
                a[k] += v
            else:
                a[k] = v
        return QuadraticRadical(*((v, k) for k, v in a.items()))

    @overload
    def __sub__(self, other: Union[
        _Integral, _Rational, float, "QuadraticRadical"
    ]) -> "QuadraticRadical":
        ...
    @overload
    def __sub__(self, other: Union[_Real, "QuadraticRadical"]) -> float:
        ...
    @overload
    def __sub__(self, other: Complex) -> complex:
        ...
    def __sub__(self, other):
        if isinstance(other, (Rational, float)):
            return self - QuadraticRadical(Fraction(*other.as_integer_ratio()))
        if not isinstance(other, QuadraticRadical):
            if isinstance(other, Real):
                return float(self) - other
            if isinstance(other, Complex):
                return complex(self) - other
            return NotImplemented
        a: Dict[int, Fraction] = self.quantity.copy()
        for k, v in other.quantity.items():
            if k in a:
                a[k] -= v
            else:
                a[k] = -v
        return QuadraticRadical(*((v, k) for k, v in a.items()))

    @overload
    def __mul__(self, other: Union[
        _Integral, _Rational, float, "QuadraticRadical"
    ]) -> "QuadraticRadical":
        ...
    @overload
    def __mul__(self, other: Union[_Real, "QuadraticRadical"]) -> float:
        ...
    @overload
    def __mul__(self, other: Complex) -> complex:
        ...
    def __mul__(self, other):
        if isinstance(other, (Rational, float)):
            return self * QuadraticRadical(Fraction(*other.as_integer_ratio()))
        if not isinstance(other, QuadraticRadical):
            if isinstance(other, Real):
                return float(self) * other
            if isinstance(other, Complex):
                return complex(self) * other
            return NotImplemented
        a: Dict[int, Fraction] = {}
        for k1, v1 in self.quantity.items():
            for k2, v2 in other.quantity.items():
                simplify_result: Tuple[int, int] = _simplify(k1 * k2)
                if simplify_result[1] in a:
                    a[simplify_result[1]] += \
                        Fraction(v1 * v2 * simplify_result[0])
                else:
                    a[simplify_result[1]] = \
                        Fraction(v1 * v2 * simplify_result[0])
        return QuadraticRadical(*((v, k) for k, v in a.items()))

    @overload
    def __truediv__(self, other: Union[
        _Integral, _Rational, float, "QuadraticRadical"
    ]) -> "QuadraticRadical":
        ...
    @overload
    def __truediv__(self, other: Union[_Real, "QuadraticRadical"]) -> float:
        ...
    @overload
    def __truediv__(self, other: Complex) -> complex:
        ...
    def __truediv__(self, other):
        if isinstance(other, (Rational, float)):
            return self / QuadraticRadical(Fraction(*other.as_integer_ratio()))
        if not isinstance(other, QuadraticRadical):
            if isinstance(other, Real):
                return float(self) / other
            if isinstance(other, Complex):
                return complex(self) / other
            return NotImplemented
        cd: QuadraticRadical = QuadraticRadical(1)
        fr: QuadraticRadical = other
        othl: List[int] = list(other.quantity)
        for i in range(1, 2**len(other.quantity)):
            if fr.rational() is not None:
                break
            r = QuadraticRadical(*((-other.quantity[othl[j]] if (i>>j)&1 else
                                    other.quantity[othl[j]], othl[j]) for j in
                                   range(len(othl))))
            fr *= r
            cd *= r
        rat: Optional[Fraction] = fr.rational()
        assert rat is not None
        return QuadraticRadical(*((v/rat, k) for k, v in
                                  self.quantity.items())) * cd

    def __floordiv__(self, other):
        a = self / other
        b = math.floor(a)
        return QuadraticRadical(b) if isinstance(a, QuadraticRadical) else b

    def __mod__(self, other):
        return self - (self // other) * other

    @overload
    def __pow__(self, power: _Integral, modulo=None) -> "QuadraticRadical":
        ...
    @overload
    def __pow__(self, power: _Real, modulo=None) -> float:
        ...
    @overload
    def __pow__(self, power: _Complex, modulo=None) -> complex:
        ...
    def __pow__(self, power, modulo=None):
        if isinstance(power, Integral):
            a: int = int(power)
            if not a:
                return QuadraticRadical(1) if modulo is None else \
                    QuadraticRadical(1) % modulo
            if a < 0:
                b = (QuadraticRadical(1)/self) ** -a
                if modulo is not None:
                    b %= modulo
                return b
            cachedbinpows: List["QuadraticRadical"] = [self]
            required: List[int] = []
            c: int = a
            for i in range(a.bit_length()-1):
                cachedbinpows.append(cachedbinpows[-1] * cachedbinpows[-1])
                if c & 1:
                    required.append(i)
                c >>= 1
            res: QuadraticRadical = cachedbinpows[-1]
            for i in reversed(required):
                res *= cachedbinpows[i]
            return res
        if isinstance(power, Real):
            return pow(float(self), power, modulo)
        if isinstance(power, Complex):
            return pow(complex(self), power, modulo)
        return NotImplemented

    def __radd__(self, other):
        return self + other

    def __rsub__(self, other):
        return -(self - other)

    def __rmul__(self, other):
        return self * other

    def __rtruediv__(self, other):
        if isinstance(other, (Rational, float)):
            return QuadraticRadical(Fraction(*other.as_integer_ratio())) / self
        if isinstance(other, Real):
            return other / float(self)
        if isinstance(other, Complex):
            return other / complex(self)
        return NotImplemented

    def __rfloordiv__(self, other):
        if isinstance(other, (Rational, float)):
            return QuadraticRadical(Fraction(*other.as_integer_ratio())) //self
        if isinstance(other, Real):
            return other // float(self)
        return NotImplemented

    def __rmod__(self, other):
        if isinstance(other, (Rational, float)):
            return QuadraticRadical(Fraction(*other.as_integer_ratio())) % self
        if isinstance(other, Real):
            return other % float(self)
        return NotImplemented

    def __rpow__(self, other):
        return other ** float(self)

    def __eq__(self, other) -> bool:
        if not isinstance(other, QuadraticRadical):
            return float(self) == other
        return self.quantity == other.quantity

    def __ne__(self, other) -> bool:
        return not self == other

    def __gt__(self, other) -> bool:
        if not isinstance(other, QuadraticRadical):
            return float(self).__gt__(other)
        return (self-other).sign() > 0

    def __ge__(self, other) -> bool:
        if not isinstance(other, QuadraticRadical):
            return float(self).__ge__(other)
        return (self-other).sign() >= 0

    def __lt__(self, other) -> bool:
        if not isinstance(other, QuadraticRadical):
            return float(self).__lt__(other)
        return (self-other).sign() < 0

    def __le__(self, other) -> bool:
        if not isinstance(other, QuadraticRadical):
            return float(self).__le__(other)
        return (self-other).sign() <= 0
