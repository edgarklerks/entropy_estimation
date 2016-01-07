{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE Unsafe  #-}
{-|
Module : Grid
Description : A grid representing a solid
Copyright : Edgar Klerks
License : BSD-3
Maintainer : Edgar.Klerks@sanoma.com
Stability : experimental
Portability : POSIX

= Introduction

This module defines a lattice of state vectors (called a 'Grid'). It represent a solid, which doesn't have interaction with its neighbouring cells.
-}
module Grid where

import qualified Data.Vector as V
import Counting


-- | = Grid
-- A 'Grid' has an associated 'StateVector', which is an abstract representation of
-- all variables of one cell.
-- All cells are contained and indexed by the 'Lattice'.
-- State vectors can be represented as Integers, so states are explicitly discrete.
class Grid f where
      -- | The lattice contains all the StateVectors, it is indexed by i, so it can be thought of as:
      -- (i, StateVector f)
      data Lattice i f :: *
      -- | A state vector contains the entire state of one cell, that is all the variables
      -- needed to define one cell.
      --
      -- == StateVector
      data StateVector f :: *
      -- | Uniquely encode a 'StateVector' as 'Integer'.
        -- The following should hold:
        --
        -- prop> a == b && encodeStateVector a == encodeStateVector b
        -- prop> a /= b && encodeStateVector a /= encodeStateVector b
        --
      encodeStateVector :: StateVector f -> Integer
      -- | A 'Lattice' is indexed by i and it indexes 'StateVectors'
      -- with 'getStateVector' and 'setStateVector' individual StateVectors could be
      -- manipulated.
      --
      -- The following should hold for 'getStateVector' and 'setStateVector':
      --
      -- prop> getStateVector (setStateVector v i l) i == l
      -- Get a individual StateVector from a Lattice
      getStateVector :: Lattice i f -> i -> StateVector f
      -- | Set a individual StateVector in a Lattice
      --
      setStateVector :: Lattice i f -> i -> StateVector f -> Lattice i f
      -- | Calculate the internal energy of a state
      --
      -- == Lattice
      --  The following methods are meant for lattices.
      --
      stateVectorEnergy :: StateVector f -> Double
      -- | Calculate the entropy of a 'Lattice'
      -- Default definitions is:
      --
      -- @
      -- latticeEntropy = uncurry entropy' . encodeLattice
      -- @
      -- Usually a user, would want to adjust the above to get the units right.
      -- See 'entropy'' to see how the entropy is calculated.
      latticeEntropy :: Lattice i f -> Double
      latticeEntropy = uncurry entropy' . encodeLattice

      -- | Encode a lattice as a whole, it will give the size of
      -- the set of all states and the states encoded as number of occurrences.
      --
      -- E.g:
      -- A lattice with alphabet "abc" and state "aaa" will become: (3, [3])
      --
      -- while "abc" "abb" will become: (3, [1,2])
      encodeLattice :: Lattice i f -> (Integer, [Integer])
