{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE Safe #-}
{-|
module : Entropy
Description : Various methods to calculate the entropy of a information theoretic solid
Copyright: Edgar Klerks
License : BSD-3
Maintainer : edgar.klerks@sanoma.com
Stability : experimental
Portability : POSIX

= Introduction

This module contains various methods of calculating the entropy of a information theoretic
solid from its state vectors. It only supports discrete states.

This is done by calculating the multiplicity of the solid from a encoded state vector.

-}
module Entropy where

-- | Count all possible permutations of a multiset, thus for aac this is:
--
-- * aac
-- * aca
-- * caa
-- * = 3
--
-- The above case should be run length encoded to: [2,1]
--
-- > multisetMultiplicity [2,1] == 3
-- > multisetMultiplicity [1,1,1] == 6
multisetMultiplicity :: [Integer]
                     -> Integer
multisetMultiplicity xs = let n = sum xs
                          in fac n `div` (product (fac <$> xs))

-- | Faculty not optimized for big numbers
--

-- | = Generalized Einstein solids
--
-- Then there is this other idea, given a set of states and cells,
-- it is possible to calculate the amount of combinations directly.
--
-- Assuming that we have states {a,b,c..} of which we have exactly
-- {an, bn, cn..} units to divide and we have n cells, we can derive an
-- expression directly. Imagine that we have 3 units of state a,
-- 5 units of state b and 4 units of state c and 3 cells. We can then
-- write a random state:
--
-- > (1,1,1) | (1,2,2) | (1,2,1)
--
-- And another one:
--
-- > (0,2,1) | (1,1,2) | (2,2,1)
--
-- Because the states are independent, we can rewrite the above as three states:
--
-- > 0 | 1 | 2
-- > 2 | 1 | 2
-- > 1 | 2 | 1
--
-- Each is an independent normal einstein solid!
--
-- Thus taking the expression for multiplicity of an einstein solid as starting point:
--
-- > (q + N - 1)! / ( q! (N - 1)!)
--
-- Where q is the number of states and N is the number of state vectors (thus N - 1 is the number of cells).
--
-- Given states {s1, s2, s3... sn} with units {q1, q2, q3 ... qn} and N cells we can state the generalized expression as follow:
--
-- > Product(i=0..n) (qi + N - 1)! / qi! (N - 1)!
--
-- Another way to view is that the multiplicitly is the same as the probability if we pick a random system from all systems. Because we assume that the states are not dependent, we can simply calculate:
--
-- > P(A,B,C) = P(A) P(B) P(C)
--
-- The next step is to get the entropy S. This is the logarithm:
--
-- > ln (Product(i=0..n) (qi + N - 1)! / qi! (N - 1)!)
--
-- Which expands to (due to log(a * b) = log a + log b
-- and log(a / b) = log a - log b
--
-- > Sum(i=0..n) ln (qi + N - 1)! - ln (qi! (N - 1)!)
--
-- Now we note that N - 1 = N, when N >> 1
--
-- > Sum(i=0..n) ln (qi + N)! - (ln (qi!) + ln (N!))
--
-- Invoking stirlings approximation:
--
-- > Sum(i=0..n) (qi + N) (ln (qi + N) - 1) - (qi ln qi - qi + N ln N - N)
-- > Sum(i=0..n) (qi + N) (ln (qi + N) - 1) - qi ( ln qi - 1) - N ( ln N - 1)
--
-- Taking the terms together (dropping the sum, because of lazy):
--
-- > qi ( ln (qi + N) - 1 - ln qi + 1) + N (ln (qi + N) - 1 - ln N + 1)
-- > qi (ln (qi + N) - ln qi) + N (ln (qi + N) - ln N)
--
-- Drawing the factors in again:
--
-- > qi (ln (qi + N)/qi) + N (ln (qi + N)/N)
-- > qi (ln (1 + N/qi) + N (ln (qi/N + 1))
--
-- Which is rather nice as result

fac :: Integer -> Integer
fac n = product [1..n]

-- | The multiplicity of a generalized einstein solid state (finally, yaj!)
--
-- = States which aren't independend
--
-- Back to the probability. If two state variables aren't independent, e.g.
-- they are describing the same underlying quantity, then genEinstein shouldn't
-- be used on those.
--
-- Say we have quantity P and Q, where P can be derived from Q, that means there is
-- a function so that:
--
-- > f :: Q -> P
--
-- Then only Q should be used to calculate the entropy.
--
-- Such a variable is called a variable of state and f is called an equation of state .
--
-- = Multiplicity
--
-- Multiplicity is the size of the state space of a system. E.g. for a 2 dimensional system with two state variables p and q, we get:
--
-- > Omega_total = Omega_p * Omega_q
--
-- Thus you only have to calculate the multiplicity for each state variable, then you can simply multiply them.
-- To get the entropy, simply dump a logarithm before it and add a constant to deal with the units (in informational solids/gasses you should end up with bits per symbol etc):
--
-- S = c log (Omega_total)
--
-- Multiplicity is not the same as the probability of a system, but it says something about how common the system is. (In Dutch it would be waarschijnlijkheid, there is no equivalent in the English language).
--
-- There is another system to be described.
--
-- Say we have a system with N particles, which are in state 1 to m, where n1 .. nm is a count of systems occupying in state 1 .. m, then we can describe it as follow:
--
-- > Omega =  N! / (Product(i=1..n) ni!)
-- > S = log (N! / (Product(i=1..n) ni!))
--
-- Using log(a/|*b) = log a -|+ log b again:
--
-- > S = log N! - Sum(i=1..n) Log (ni!)
--
-- Invoking stirling:
--
-- > S = N log N - N - Sum (ni log ni - ni)
--
genEinstein :: Double -- ^ Number of cells
            -> [Double] -- ^ Number of units available
            -> Double -- ^ Entropy
genEinstein n qi  = sum (fmap calc qi)
            where calc q = q * log (1 + n/q) + n * log (1 + q/n)

-- | Entropy for a system of n particles, which can be in s1 .. sm.
-- genEinstein should be used if you know the total quantity of the state variables, e.g 30 units of energy over 30 cells.
-- genMultiset should be used if you know how much particles are in each state, e.g. 20 are in state a, 30 are in state b, 40 are in state c etc.
genMultiSet :: Double -- ^ Number of state vectors
            -> [Double] -- ^ Number of particles in a particular state.
            -> Double
genMultiSet n qs = n * log n - n - sum ((\q -> q * log q - q) <$> qs)
-- | Stirling approximation of faculty, accurate for big numbers.
--
-- /Note/ that we are dropping into floating territory.
facStir :: Double -> Double
facStir p = sqrt (2 * pi * ee * p) * (p/ee)**p

-- | Stirling approximation for the factorial function under a logarithm
-- Becomes more and more accurate over p>100
logfacStir :: Double -> Double
logfacStir p = p * log p - p



-- | The value /e/
ee :: Double
ee = exp 1
