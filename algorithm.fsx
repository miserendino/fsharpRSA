open System
open System.Numerics

// In F# sqrt non è implementa per il tipo BigInteger, l'implementazione deriva da http://www.codeproject.com/Articles/266425/BigInteger-Square-Root-in-Fsharp 

//*******************************************************************************************************************

// Here are some constant BigIntegers to help with the log10 function.

let ten = 10I  // 10^1

let tenBillion = 10000000000I // 10^10

let googol = 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000I // 10^100

let googolToTenth = googol * googol * googol * googol * googol * googol * googol * googol * googol * googol // 10^10000



//*************************************************************************************

// Take the integer logarithm of theNumber, base 10. 

let log10 theNumber =

    if theNumber <= 0I then

        failwith "Logarithm of a non-positive number is not defined."

    // this inner functions counts the number of times divisor will divide num

    let rec divisions count (num : BigInteger) divisor =

        match num with          // stiamo semplicemente vedendo quante volte il divisor sta nel num e incrementiamo count ogni volta

        | x when x < divisor -> count, num

        | _ -> divisions (count + 1I) (num / divisor) divisor

    let thousandsCount, rem1 = divisions 0I theNumber googolToTenth

    let hundredsCount, rem2 = divisions 0I rem1 googol

    let tensCount, rem3 = divisions 0I rem2 tenBillion

    1000I * thousandsCount + 100I * hundredsCount + ten * tensCount + fst(divisions 0I rem3 10I)





//*************************************************************************************

// Raises 10I to the power of a specified value.

// Warning: Be careful with passing in a value over 100000I (hundred thousand).

//          It could take a long time to complete.

// (BigInteger -> BigInteger)

let private power10 (exponent : BigInteger) =

    if exponent.Sign = -1 then

        0I

    else

        let rec expox (product : BigInteger) (pow : BigInteger) multiplier log10OfMult =

            match pow with

            | x when x < log10OfMult -> product, pow

            | _ -> expox (product * multiplier) (pow - log10OfMult) multiplier log10OfMult

        let pow10To1000, rem1 = expox 1I exponent googolToTenth 1000I

        let pow10To100, rem2 = expox 1I rem1 googol 100I

        let pow10To10, rem3 = expox 1I rem2 tenBillion 10I

        pow10To1000 * pow10To100 * pow10To10 * fst(expox 1I rem3 10I 1I)





//*************************************************************************************

// Take the integer square root of a specified value.

// Because this is dealing with integers the answere is truncated.

// An exception occurs if you pass in a negative value.

// (BigInteger -> BigInteger)

let bigintSqrt(bigNum : BigInteger) =

    let sqrtRoughGuess (num : BigInteger) =                 

        let log10x = log10 (num + 1I)

        let halfLog = (log10x + 1I) >>> 1       // This is because 10 is 2 in binary. Multiplying a number by 10(be it binary or decimal or hexadecimal) appends a 0 to the number(which is effectively left shifting). Similarly, dividing by 10(or 2) removes a binary digit from the number(effectively right shifting). 

        (power10 halfLog)

    let rec converge prevGuess =

        let nextGuess = (bigNum / prevGuess + prevGuess) >>> 1      // 1/2 ( S/x + x )

        match BigInteger.Abs (prevGuess - nextGuess) with

        | x when x < 2I -> nextGuess            // precisione arbitraria

        | _ -> converge nextGuess

    if bigNum.Sign = -1 then

        failwith "Square root of a negative number is not defined."

    else

        let root = converge (sqrtRoughGuess bigNum)     // prendiamo come valore inziale per l'agoritmo il logaritmo in base 10 del numero + 1 diviso 2 elevato alla 10

        if root * root > bigNum then

            root - 1I

        else

            root

// this verifies if input is a perfect root

let IsPerfectSquare (input : BigInteger) = 
  let closestRoot = bigintSqrt input
  input = closestRoot * closestRoot;

// same as mpz_sqrtrem, calcola la radice dell'input poi capisce se è una radice perfetta

let sqrtrem (input : BigInteger) = 
  let rop1 = bigintSqrt input
  let rop2 = input - rop1 * rop1
  rop2
           
      
//*************************************************************************************************

// se è una radice perfetta allora restituiscimi a

// con accumulatore 

let rec loopF (a : BigInteger) (v : BigInteger) = 
        match (sqrtrem (a * a - v)) with 
          | x when x = 0I -> a 
          | _ -> loopF (a + 1I) v

// Fermat's factorization method 

let fermat (n : BigInteger)= 
  let a = if ( (sqrtrem n) = 0I ) then bigintSqrt n else (bigintSqrt n) + 1I  // ceiling (sqrt N), parte intera superiore
  let aperamenon = BigInteger.Pow(loopF a n, 2 ) - n        // a^2 - n
  let b = bigintSqrt aperamenon     // b = sqrt(b2)
  printf "p = %A\nq =  %A\n" (a + b) (a - b)    // il numero n è composta da (a + b)*(a - b)

let _ = fermat 7723537042671053533349250770213334979939936820556134073043654409875795653759249636472530580914872873851196626721386328853913207026848732669575096806570611I;;

// find the private exponent d (multiplicative inverse in mod phi(n) )

exception OutOfBoundException

let privExp (p : BigInteger) (q : BigInteger) (e : BigInteger) = 
  let phi_n = (p - 1I) * (q - 1I) in 
    let rec loop k =
        match BigInteger.Remainder(phi_n * k , e) with      // trova i valori per cui d = 1 + k(p-1)(q-1)
          | x when x = 0I -> (phi_n * (k - 1I) ) / e
          | x when x < e -> loop (k + 1I)  
          | _ -> raise OutOfBoundException      
    loop 2I
    
let _ = privExp 148014008892085729643675437352771337123457425613021769133485669342043764166248475537091040732725971887614029502547017118261383577275568664346427549123275573899584262102353725978410291261322011923743903081324644222331364445917991336631407165072327107290181592383959906391578626558349562465943656697778534442483I 148014008892085729643675437352771337123457425613021769133485669342043764166248475537091040732725971887614029502547017118261383577275568664346427549123275573899584262102353725978410291261322011923743903081324644222331364445917991336631407165072327107290181592383959906903120796472368462831702487846238746829211I 3I;;

// python implementation

(*
def privateExponent(p,q,e):
    phi_n=(p-1)*(q-1)
    for k in range(1,e):            
        if (phi_n*k+1) % e==0:
            return (phi_n*k+1)/e
    return -1 
*)

// better with continuations

let id x = x;

let factorialB (n : BigInteger) = 
  let rec loop (n : BigInteger) k = 
    match n with 
      | n when n = 1I -> k 1I
      | n when n > 0I -> loop (n - 1I) (fun res -> k ( res * n))
      | _ -> failwith "Negative number factorial"
  loop n id


// Pollard p-1 algorithm

let pollard a limit n =
  let rec loop (acc1 : BigInteger) acc2 = 
    match acc2 with 
      x when x = limit -> BigInteger.Remainder(factorialB a , n) 
      | x when x < limit -> if ( BigInteger.GreatestCommonDivisor(acc1 - 1I,n) > 1I && BigInteger.GreatestCommonDivisor(acc1 - 1I,n) < n ) 
                            then BigInteger.GreatestCommonDivisor(acc1 - 1I,n)              // se MCD(b-1,n) abbiamo trovato un fattore non banale di n
                            else loop ( BigInteger.Remainder( BigInteger.Pow(acc1 , acc2) , n ) ) (acc2+1)          // il primo acc contiene le potenze crescenti il secondo l'accumulatore
      | _ -> failwith "An error occured" 
  loop ( BigInteger.Remainder(a,n) ) 1          // iniziamo con b1 = a (mod n)

let _ = pollard 2I 100000 51069631I;;

let _ = pollard 2I 100000 (100000004021I * 100000006409I * 100000007837I);;

// Miller-Rabin primality test

let scompose n =  
  let m = n - 1I in
    let rec loop t s = 
        match (t % 2I) with            
            | x when x = 0I -> loop (t >>> 1) (s + 1I)
            | x when x <> 0I -> t, s                // qualunque numero dispari lo posso scrivere come 2^s * t + 1
            | _ -> failwith "input error"    
    loop m 0I

let _ = scompose 189I;;
let _ = scompose 2349872034650234650136509134650246026450872346501365I;;

let BigInt_pow (bas : BigInteger) (exp : BigInteger) =
    let rec loop acc1 acc2  = 
      match acc1 with 
          | x when x = 0I -> acc2
          | x when x > 0I -> loop (acc1 - 1I) (acc2 * bas)
          | _ -> failwith "Error, negative exponent"
    loop exp 1I

let inner_milRab b r t n = 
    let rec loop acc1 acc2 =
      match (acc1 % n) with 
        | x when x = -1I -> printf "%A è un primo o pseudoprimo forte su base %A\n" n b
        | _ -> match acc2 with 
                             | y when y < r -> loop (BigInt_pow b ( (BigInt_pow 2I acc2) * t) ) (acc2 + 1I)
                             | _ -> printf "%A non è primo\n" n
    loop (BigInt_pow b ( ( BigInt_pow 2I 1I ) * t) ) 1I

let miller_rabin (a : BigInteger) (n : BigInteger) = 
  let m,k = scompose n in let b = (BigInt_pow a m) % n   
      if (b = 1I) then
                     let b1 = (BigInt_pow b m) % n 
                     if ( b1 = 1I || b1 = -1I) then printf "%A è un primo o pseudoprimo forte su base %A\n" n b
                                                else 
                                                 let r = k - 1I 
                                                 inner_milRab b r m n
                  else
                  printf "%A non è primo\n" n

// Solovay-Strassen algorithm

let rec gcd a b =
      if b = 0I then a
      else gcd b (a % b)

// jacobi symbol is defined for all odd integers a

let rec jacobi a n =
    if (a = 0I && (n = 1I)) then 1I elif ( a = 0I && (n <> 1I) ) then 0I
    elif ( a = 2I ) then 
      match ( n % 8I ) with 
        | x when x = 1I || x = 7I -> 1I
        | x when x = 3I || x = 5I -> -1I
        | _ -> failwith "Error"
    elif ( a >= n ) then jacobi (a % n) n
    elif ( a % 2I = 0I ) then (jacobi 2I n) * (jacobi (a >>> 1) n)
    elif ( a  % 4I = 3I && n % 4I = 3I ) then -(jacobi n a) 
    else jacobi n a
    
let _ = jacobi 6278I 9975I;;        // -1

let _ = jacobi 1236I 20003I;;       // 1 

// a deve essere minore di n e i due devono essere coprimi tra loro  

let rec goodAChosen n = 
  let rnd = System.Random()
  let a = BigInteger ( rnd.Next() ) in 
    match gcd a n with 
        | x when x = 1I -> a
        | _ -> goodAChosen n


let sol_str n (iter : int) =
  let rec loop acc = 
    match acc with 
        | x when x < iter -> let a = goodAChosen n in            // genera int32  
                             let j = ( n + jacobi a n ) % n
                             if ( j <> 0I || ( BigInt_pow a (( n - 1I) >>> 1) ) % n <> j) 
                             then loop (acc + 1)
                             else printf "%A è probabilemente primo\n" n
        | _ -> printf "%A è probabilemente primo\n" n
  loop 0
               
// questo test presenta una precisione imposta arbitrariamente affinchè l'algoritmo sia più accurato

let solovay_strassen n =
    if ( n < 2I ) then printf "%A non è primo\n" n
    elif ( n % 2I = 0I && n <> 2I ) then printf "%A non è primo\n" n
    else 
    sol_str n 50     // precisione imposta arbitrariamente





