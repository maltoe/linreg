type Vector = [Double]
type Matrix = [[Double]]

dim :: Matrix -> (Int, Int)
dim [] = (0, 0)
dim [x] = (1, length x)
dim (x:xs)
    | length x /= c = error "dim: Corrupt matrix."
    | otherwise = (r + 1, c)
    where
    (r, c) = dim xs

transpose :: Matrix -> Matrix
transpose ([]:_) = []
transpose x = ((map head x):(transpose (map tail x)))

multiply :: Matrix -> Matrix -> Matrix
multiply a b
    | acols /= brows = error "multiply: Invalid dimensions."
    | otherwise = [[ sum (zipWith (*) av bv) | bv <- bt ] | av <- a]
    where
    (_, acols) = dim a
    (brows, _) = dim b
    bt = transpose b
    
add :: Matrix -> Matrix -> Matrix
add [] [] = []
add (a:as) (b:bs)
    | dim (a:as) /= dim (b:bs) = error "add: Invalid dimensions."
    | otherwise = (zipWith (+) a b):(add as bs)
    
scale :: Double -> Matrix -> Matrix
scale _ [] = []
scale l (a:as) = ((map (* l) a):(scale l as))
       
det :: Matrix -> Double
det [] = error "det: Empty matrix."
det [[a]] = a
det a
    | (rows /= cols) = error "det: Invalid dimensions."
    | otherwise = sum [ (-1)^(j + 1) * (head a)!!(j - 1) * det (_remove a 0 (j - 1)) | j <- [1..cols] ]
    where
    (rows, cols) = dim a
    
cofactorMatrix :: Matrix -> Matrix
cofactorMatrix a
    | cols /= rows = error "cofactorMatrix: Invalid dimensions."
    | otherwise = [[ cofactor a i j | j <- [1..cols] ] | i <- [1..rows]]
    where
    (rows, cols) = dim a
    cofactor x y z = (-1)^(y + z) * det (_remove x (y - 1) (z - 1))
    
inverse :: Matrix -> Matrix
inverse a = transpose [[ x / (det a) | x <- row ] | row <- (cofactorMatrix a) ]

_remove :: Matrix -> Int -> Int -> Matrix
_remove a i j = transpose (removeRow (transpose (removeRow a i)) j)
    where
    removeRow x n = take n x ++ drop (n + 1) x

k = [[1.0,2.0],[3.0,4.0]]

polynomialFeatures :: Int -> Matrix -> Matrix
polynomialFeatures n a = _prependOne (map (\av -> foldl1 (++) [ map (foldl1 (*)) (_kcombination av i) | i <- [1..n]]) a) 

linearFeatures = polynomialFeatures 1
quadraticFeatures = polynomialFeatures 2
cubicFeatures = polynomialFeatures 3

_prependOne :: Matrix -> Matrix
_prependOne [] = []
_prependOne (a:as) = ((1.0:a):(_prependOne as))

_kcombination :: [a] -> Int -> [[a]]
_kcombination _ 0 = [[]]
_kcombination [] _ = []
_kcombination (a:as) k = map (a:) (_kcombination (a:as) (k - 1)) ++ _kcombination as k

solveForBeta :: Matrix -> Vector -> (Matrix -> Matrix) -> Matrix
solveForBeta x y phi = multiply (multiply (inverse (multiply phixt phix)) phixt) yt
    where
    phix = phi x
    phixt = transpose phix
    yt = transpose [y]
    
_parseFile :: String -> String -> [Double] -> [[Double]]
_parseFile [] part acc
    | length part /= 0 = [(acc ++ [read part :: Double])]
    | length part == 0 = [acc]
_parseFile (x:xs) part acc
    | (x == ' ') && ((length part) /= 0) = _parseFile xs "" (acc ++ [read (part ++ [x]) :: Double])
    | (x == ' ') && ((length part) == 0) = _parseFile xs part acc
    | (x == '\n') && ((length part) /= 0) = [acc ++ [read (part ++ [x]) :: Double]] ++ _parseFile xs "" []
    | (x == '\n') && ((length part) == 0) = [acc] ++ _parseFile xs part []
    | otherwise = _parseFile xs (part ++ [x]) acc

_extractX :: [[Double]] -> Matrix
_extractX [] = []
_extractX (x:xs) = ((take ((length x) - 1) x):(_extractX xs))

_extractY :: [[Double]] -> Vector
_extractY [] = []
_extractY (x:xs) = (drop ((length x) - 1) x) ++ _extractY xs

main = do
    s <- readFile "dataLinReg1D.txt"
    print (_extractY (_parseFile s "" []))
    print (solveForBeta (_extractX (_parseFile s "" [])) (_extractY (_parseFile s "" [])) linearFeatures)
    
