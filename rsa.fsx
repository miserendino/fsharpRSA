open System
open System.Numerics

let rec BigInt_pow (bas : BigInteger) (exp : BigInteger) = 
    match exp with 
        | x when x = 0I -> 1I
        | x when x > 0I -> bas * BigInt_pow bas (exp - 1I)
        | _ -> failwith "Error, negative exponent"


let rsa p q m e = 
  let n = p * q
  (BigInt_pow m e) % n

let rec gcd a b = 
  if b = 0I then a else gcd b (a % b) 

// using extended Euclid algorithm

let mod_inv a n =
       let rec ext_euclid a b =
          if (b = 0I) then (a, 1I, 0I) 
          else let (d', x', y') = ext_euclid b (a % b)  // prendiamo solo il resto della divisione intera
               let (d, x, y) = (d', y', x' - (a/b)*y')  
               (d, x, y)

       let (d, x, y) = ext_euclid a n
       if (x < 0I) then (x+n) else x

let decrypt p q c e = 
  let n = p * q
  let phi_n = (p - 1I) * (q - 1I)
  let d = mod_inv e phi_n
  (BigInt_pow c d) % n


// naive example 

let c = rsa 11I 3I 4I 3I        // 31

let d = decrypt 11I 3I c 3I     // 4


    