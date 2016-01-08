{-# LANGUAGE GADTs #-}
{-# LANGUAGE ViewPatterns #-}
{- *
 This module contains an entropy estimator, which uses statistical entropy encoding for the estimation.

 I haven't looked into dictionary type encodings, because they are not treating strings not equally. E.g.
 A source produces:

   ababababababab
   aaabaabaaaaabb

Then the dictionary generator would pick up ab as symbol in the first case and
some wild thing in the second case. This give incomparable estimates.

I am interested in trying a grammar based method, which tries to search for the smallest grammar
For this we need to find a straight line grammar (which is a start symbol s and some production rules P1..Pn, which
generates exactly 1 string).


-}
module Main where

import Data.List
import Debug.Trace
import qualified Data.Vector as V
import Data.Maybe
import qualified Data.Map.Strict as S
import Data.Ratio
import Data.Word
import Control.Applicative

type Model a = a -> (Rational, Rational)

type Prob a = [(a, Rational)]

type StaticModel a = (S.Map a Rational, S.Map a Rational)

data EncoderState a = EC {
                  currentCodePoint :: Rational,
                  currentIntervalWidth :: Rational,
                  getProb :: Model a
        }

instance Show (EncoderState t) where
  show (EC lr hr _) = "Interval " ++ (show (lr, hr + lr))

{-- Entropy encoder --}
encodedInterval :: EncoderState t -> (Rational, Rational)
encodedInterval (EC lr hr _) = let rr = lr in (rr, hr + rr)

startState :: Model a -> EncoderState a
startState  = EC 0 1

encode :: Ord a => [a] -> Model a -> EncoderState a
encode xs m = foldl' step ( startState m) ( xs)
       where step z x = let iwidth = currentIntervalWidth z
                            icpoint = currentCodePoint z
                            (prob, cprob) = getProb z x
                        in z {
                           currentIntervalWidth = prob * iwidth,
                           currentCodePoint = icpoint + cprob * iwidth
                        }
-- | Optimized bwt transform
-- Space O(n), Time O(n)
-- The space efficiency is achieved by cycle the list, and calculate the position of the last column and store it together with the first column,
-- then sort it and shave of the first column. Time is actually N^2 because we are using lists, but with arrays it would be N
bwt :: Ord a => [a] -> [a]
bwt xs = let ys = Nothing : Just `fmap` xs
         in let lx = length ys
            in strip $ sortBy sortf $ worker lx lx (cycle ys)
    where worker jmp 0 _ = []
          worker jmp n (x:xs) = let ns = drop jmp xs
                                      in (x, head $ ns ) : worker jmp (n - 1) ns
          sortf (fst -> a) (fst -> b) | a == Nothing = LT
                                      | b ==  Nothing = GT
                                      | otherwise = compare a b
          strip [] = []
          strip (x:xs) = let (a,b) = x in case b of
                                            Nothing -> strip xs
                                            Just s -> s : strip xs

-- | Move to front transform
mtf :: Ord a => [a] -> [a] -> [a]
mtf prob xs = let ss = fmap head (group $ sort prob)
              in toSym ss $ fst $ foldl' step ([],ss) xs
      -- index is partial
      where index a ss = (fromJust $ findIndex (==a) ss, takeWhile (/=a) ss, dropWhile (/=a) ss)
            step (z,ss) x = let (i, xs,y:ys) = index x ss
                              in  (i:z, y: (xs ++ ys) )

toSym :: Ord a => [a] -> [Int] -> [a]
toSym prob xs = fmap (prob !!) xs

{- Model building functions --}
staticModel :: Ord a => StaticModel a -> Model a
staticModel m a =  ( ( fst m) S.! a,  (snd m )S.! a)

idealModel :: Ord a => [a] -> Model a
idealModel  = probToModel . toProb . sortBy sortf . fmap getL . group . sort
           where getL xs = (head xs, fromIntegral $ length xs)
                 sortf a b = snd b `compare` snd a
                 toProb xs = let l = sum (snd <$> xs)
                             in (\(a,p) -> (a, p / l :: Rational) ) <$> xs

probToModel :: Ord a => Prob a -> Model a
probToModel xs = staticModel $ cumProbs xs

cumProbs :: Ord a => Prob a -> StaticModel a
cumProbs xs = let p = fst $ foldlKeys step (S.empty, 0) xs in (S.fromList xs, p)
              where step k v (sm, c) = let c' = c + v
                                       in (S.insert k c sm, c')

foldlKeys :: (a -> b -> c -> c) -> c -> [(a,b)] -> c
foldlKeys step z xs = foldl' step' z xs
          where step' z (a,b) = step a b z

testModel :: Prob Char
testModel =  [('a', 1/2),('b', 1/4), ('c', 1/8), ('d', 1/8)]

testSequence :: String
testSequence = "abaabcda"


between :: Double -> (Double, Double) -> Bool
between x pp = fst pp <= x && x < snd pp

{-- Entropy estimation functions --}
relevantBits :: (Rational, Rational) -> Integer
relevantBits (l,r) = worker l r 0
  where worker l r n | headI l == headI r = worker (tailI l) (tailI r) ( 10 * n + headI l)
                     | headI l < headI r = 10 * n + ( headI l + 1)
                     | otherwise = 10 * n + (headI l + 1)


-- | Expects x to be between 0 <= x < 1
headI :: Rational -> Integer
headI x = truncate (x * 10)

tailI :: Rational -> Rational
tailI x = ( x - fromInteger (headI x) / 10 ) * 10

integerRep :: Ord a => [a] -> Model a -> Integer
integerRep xs m = let rr = encodedInterval $ encode xs m
                  in relevantBits rr

bitsInInteger :: Integer -> Double
bitsInInteger n = log (fromInteger n) / log 2

entropyFromModel :: Ord a => ([a] -> [a]) -> [a] -> Model a -> Double
entropyFromModel pre xs m = let p = bitsInInteger $ integerRep (pre xs) m
                            in p / (fromIntegral $  length xs)

-- | Estimate the zeroth order entropy, this is somewhat accurate if the source is long enough
idealEntropy :: (Ord b,  Ord a) => ([a] -> [b]) -> [a] -> Double
idealEntropy pre xs = let bb = pre xs in entropyFromModel id bb (idealModel bb)
