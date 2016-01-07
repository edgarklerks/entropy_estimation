{-# LANGUAGE Safe #-}
{-# LANGUAGE ViewPatterns #-}
{-|
Module : Counting
Description : Several methods to determine the entropy of a datasource
Copyright : Edgar Klerks
License : BSD-3
Maintainer : Edgar.Klerks@sanoma.com
Stability : experimental
Portability : POSIX

This module contains a method for estimating the naked entropy of a multiset. It uses the
multinomal method to count the microstates, which belong to a macrostate.

-}

module Counting where


import qualified Data.Set as S
import qualified Data.List as L

-- | Count all permutations without repetitions of a microstate
-- Output will always be a integer, but because of the division,
-- and later usage it is nicer to represent them as 'Rational'.
countMultiSetPermutations :: Ord a =>
                             [a] -- ^ Microstate
                          -> Rational
countMultiSetPermutations xs = let ss = L.sort xs
                                   ms = L.genericLength <$> L.group ss
                                   n = L.genericLength xs
                               in fac n / (product (fac <$> ms))
                          where fac n = product [1..n]


-- | Given a macro state of a system, calculate the probability
-- It is assumed that position doesn't affect the entropy, thus:
--  aba and aab are belonging to the same macrostate
probability :: (Enum a, Ord a) =>
               (a, a) -- ^ Range of states
            -> [a] -- ^ Microstate belonging to the macro state
            -> Rational -- ^ Probability
probability (l, h) rs = let ln = fromIntegral $ fromEnum h - fromEnum l + 1
                            n = L.genericLength rs
                            q = ln^n
                         in toRational (countMultiSetPermutations rs) / toRational q

-- | This will calculate the entropy when position in the list is not part of the state vector, thus
-- aba and aab refer to the same macro state.
-- This gives back the naked entropy. It could be a user needs to multiply it by a constant
-- depending on the choosen units.
-- The base for the logarithm is e.
entropy :: (Enum a, Ord a) => (a, a) -> [a] -> Double
entropy range seq = let p = fromRational $ probability range seq in - p * ( log p)

-- | A variant of 'probability', where only the size of the range of the state vector is given
-- And the state vector is encoded as number of occurences
probability' ::
                Integer -- ^ Number of states
             -> [Integer] -- ^ Occurences of states in the vector
             -> Rational -- ^ Probability
probability' ln xs = let n = sum xs
                         q = toRational $ ln^n
                         m = (toRational $ fac n) / toRational (product (fac  <$> xs))
                     in m / q
             where fac n = product [1..n]
-- | A variant of 'entropy', where only the size of the range of the state vector is given
-- and the state vector itself is encoded as number of occurences
entropy'  ::
             Integer -- ^ Number of states
          -> [Integer] -- ^ Occurences of states in the vector
          -> Double -- ^ Calculated entropy
entropy' range seq = let p = fromRational $ probability' range seq in - p * (log p )

-- | Entropy calculated using stirlings method, use when there are a lot of particles
stirlingEntropy ::
                   Integer -- ^ Number of states
                -> [Integer] -- ^ Occurences of states in the vector
                -> Double
stirlingEntropy range seq = - logprob * fromRational prob
                -- log (fac n / product (fac <$> ms) * ln^n)
                -- log (fac n) - log (fac n1 * .. fac nn) - log(ln^n)
                -- log (fac n) - (log (fac n1) + log (fac n2) .. log (fac nn)) - log(ln^n)
                where n' = fromInteger $ n
                      q' = fromInteger q
                      n = sum seq
                      q = range ^ n
                      logfac n = n * log n - n
                      fac n@(toRational -> n') = toRational (cs n') * (n' / expn) ^ n
                      cs n = dsqrt (2 * pin * expn *  n)
                      pin = toRational pi
                      expn = toRational (exp 1)
                      dsqrt n = round ( sqrt $ fromRational n)
                      prob =  ( fac $   n) / (product (fac <$> seq)) /  toRational q
                      logprob = logfac  n' - ( sum $ logfac . fromInteger <$> seq) - log q'
