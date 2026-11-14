/-
  Finite and arithmetic statements used in the paper, checked by the Lean 4 kernel.
  Core Lean only (no Mathlib). Build: lean Sandpile.lean
-/

set_option maxRecDepth 1000000
set_option maxHeartbeats 0

/-! ## The function G of Theorem B -/

/-- product of floor(l/2) and of floor((l-1)/2) over the sides -/
def P : List Nat → Nat
  | [] => 1
  | l :: ls => l / 2 * P ls

def Q : List Nat → Nat
  | [] => 1
  | l :: ls => (l - 1) / 2 * Q ls

def G (ls : List Nat) : Nat := P ls + Q ls

/-- Lemma Grec: G obeys the recursion of the layer construction. -/
theorem G_layer (l : Nat) (ls : List Nat) (h : 1 ≤ l) :
    G ((l + 2) :: ls) = G (l :: ls) + G ls := by
  have h1 : (l + 2) / 2 = l / 2 + 1 := by omega
  have h2 : (l + 2 - 1) / 2 = (l - 1) / 2 + 1 := by omega
  simp only [G, P, Q, h1, h2, Nat.add_mul, Nat.one_mul]
  omega

theorem G_one (ls : List Nat) : G (1 :: ls) = 0 := by
  simp [G, P, Q]

theorem G_two (ls : List Nat) : G (2 :: ls) = P ls := by
  simp [G, P, Q]

/-- the two-dimensional case is Dhar's formula ceil((l-1)(m-1)/2) -/
def dharOK (n : Nat) : Bool :=
  (List.range n).all fun l => (List.range n).all fun m =>
    G [l + 1, m + 1] == (l * m + 1) / 2

theorem G_dhar_upto_40 : dharOK 40 = true := by decide +kernel

theorem G_values :
    G [2, 4, 4] = 4 ∧ G [4, 4, 4] = 9 ∧ G [4, 4, 5] = 10 ∧ G [5, 5, 5] = 16 ∧
    G [6, 6, 6] = 35 ∧ G [4, 4, 4, 4] = 17 ∧ G [9, 8] = 28 := by decide

/-! ## The same-parity family of cells -/

def tuples : List Nat → List (List Nat)
  | [] => [[]]
  | n :: ns => (List.range n).foldr (fun i acc => (tuples ns).map (fun t => i :: t) ++ acc) []

/-- base points of the unit cells of the box with sides ls (0-indexed) -/
def cells (ls : List Nat) : List (List Nat) := tuples (ls.map (· - 1))

def sameParity (c : List Nat) : Bool :=
  c.all (fun x => x % 2 == 0) || c.all (fun x => x % 2 == 1)

def family (ls : List Nat) : List (List Nat) := (cells ls).filter sameParity

/-- number of common sites of the cells with base points c and c' -/
def common : List Nat → List Nat → Nat
  | x :: xs, y :: ys =>
      (if x == y then 2 else if x + 1 == y || y + 1 == x then 1 else 0) * common xs ys
  | _, _ => 1

def pairsAtMost (k : Nat) : List (List Nat) → Bool
  | [] => true
  | c :: cs => cs.all (fun c' => common c c' ≤ k) && pairsAtMost k cs

def familyOK (ls : List Nat) : Bool :=
  (family ls).length == G ls && pairsAtMost 1 (family ls)

def sortedBoxes : Nat → Nat → Nat → List (List Nat)
  | 0, _, _ => [[]]
  | d + 1, lo, hi =>
      ((List.range (hi + 1)).filter (lo ≤ ·)).foldr
        (fun l acc => (sortedBoxes d l hi).map (fun t => l :: t) ++ acc) []

/-- Lemma parity and Theorem family on all boxes with d = 2 and sides ≤ 7, d = 3 and sides ≤ 5, d = 4 and sides ≤ 3 -/
theorem family_d2 : (sortedBoxes 2 1 7).all familyOK = true := by decide +kernel
theorem family_d3 : (sortedBoxes 3 1 5).all familyOK = true := by decide +kernel
theorem family_d4 : (sortedBoxes 4 1 3).all familyOK = true := by decide +kernel

/-- the four corner cells of [2]×[4]×[4] used in the footnote are pairwise disjoint -/
theorem corner_cells_disjoint :
    pairsAtMost 0 [[0, 0, 0], [0, 0, 2], [0, 2, 0], [0, 2, 2]] = true := by decide +kernel

/-! ## Legality of the elementary histories of Section upper -/

def setAt : List Int → Nat → Int → List Int
  | [], _, _ => []
  | _ :: xs, 0, v => v :: xs
  | x :: xs, n + 1, v => x :: setAt xs n v

def getAt (xs : List Int) (n : Nat) : Int := xs.getD n 0

/-- add a grain at site s -/
def addG (h : List Int) (s : Nat) : List Int := setAt h s (getAt h s + 1)

/-- topple s, sending grains to the listed neighbours inside the box; none if illegal -/
def toppleTo (h : List Int) (s : Nat) (nbrs : List Nat) : Option (List Int) :=
  if getAt h s ≥ 2 then
    some (nbrs.foldl (fun h' t => addG h' t) (setAt h s (getAt h s - 2)))
  else none

def gray (i : Nat) : Nat := i ^^^ (i >>> 1)

/-- Lemma cube: the Gray-cycle history on [2]^k, started from the full configuration
    plus a grain at the origin, is legal and ends with a single grain. -/
def grayRun (k : Nat) : Option (List Int) :=
  let n := 2 ^ k
  (List.range n).foldl
    (fun acc j => acc.bind fun h =>
      toppleTo h (gray j) [gray j ^^^ (gray j ^^^ gray ((j + 1) % n))])
    (some (addG (List.replicate n 1) 0))

def grayOK (k : Nat) : Bool :=
  match grayRun k with
  | some h => h.foldl (· + ·) 0 == 1 && h.all (fun x => 0 ≤ x && x ≤ 1)
  | none => false

theorem gray_cube : (List.range 7).all (fun k => grayOK (k + 1)) = true := by decide +kernel

/-- Lemma chain: on [l], add a grain at the end and topple from the end to the start. -/
def chainRun (l : Nat) : Option (List Int) :=
  (List.range l).foldl
    (fun acc i => acc.bind fun h =>
      let s := l - 1 - i
      toppleTo h s ((if s > 0 then [s - 1] else []) ++ (if s + 1 < l then [s + 1] else [])))
    (some (addG (List.replicate l 1) (l - 1)))

def chainOK (l : Nat) : Bool :=
  match chainRun l with
  | some h => h.foldl (· + ·) 0 == (l : Int) - 1 && h.all (fun x => 0 ≤ x && x ≤ 1)
  | none => false

theorem chain_upto_30 : (List.range 30).all (fun i => chainOK (i + 2)) = true := by decide +kernel

/-! ## The Oslo model -/

/-- the matrix D of the path with a reflecting end, 1-indexed -/
def Dm (L x y : Nat) : Int :=
  if x == y then (if x < L then 2 else 1)
  else if x + 1 == y || y + 1 == x then -1 else 0

def DinvOK (L : Nat) : Bool :=
  (List.range L).all fun i => (List.range L).all fun j =>
    let x := i + 1; let y := j + 1
    ((List.range L).foldl (fun s k => s + Dm L x (k + 1) * (Nat.min (k + 1) y : Int)) 0)
      == (if x == y then 1 else 0)

/-- equation Dinv: D^{-1} = (min(x,y)) -/
theorem Dinv_upto_12 : (List.range 12).all (fun i => DinvOK (i + 1)) = true := by decide +kernel

def weightOK (L : Nat) : Bool :=
  (List.range L).all fun j =>
    ((List.range L).foldl (fun s i => s + ((i + 1 : Nat) : Int) * Dm L (i + 1) (j + 1)) 0)
      == (if j + 1 == L then 1 else 0)

/-- equation weightD: sum_x x D_{xy} = [y = L] -/
theorem weight_upto_12 : (List.range 12).all (fun i => weightOK (i + 1)) = true := by decide +kernel

/-- membership in R: every 0 is followed, after a run of 1's, by a 2 -/
def okAfter : List Nat → Bool
  | [] => false
  | 1 :: xs => okAfter xs
  | 2 :: _ => true
  | _ => false

def inR : List Nat → Bool
  | [] => true
  | 0 :: xs => okAfter xs && inR xs
  | _ :: xs => inR xs

def stable : Nat → List (List Nat)
  | 0 => [[]]
  | n + 1 => (stable n).foldr (fun t acc => (0 :: t) :: (1 :: t) :: (2 :: t) :: acc) []

def fib : Nat → Nat
  | 0 => 0
  | 1 => 1
  | n + 2 => fib n + fib (n + 1)

/-- Proposition count, checked by enumeration for L ≤ 7 -/
theorem countR_upto_7 :
    (List.range 7).all (fun i => ((stable (i + 1)).filter inR).length == fib (2 * (i + 1) + 1))
      = true := by decide +kernel

/-- the automaton of the proof of Proposition count -/
def ab : Nat → Nat × Nat
  | 0 => (0, 1)
  | n + 1 => (2 * (ab n).1 + (ab n).2, (ab n).1 + (ab n).2)

theorem ab_fib (n : Nat) : (ab n).1 = fib (2 * n) ∧ (ab n).2 + (ab n).1 = fib (2 * n + 1) := by
  induction n with
  | zero => exact ⟨rfl, rfl⟩
  | succ n ih =>
    obtain ⟨h1, h2⟩ := ih
    have e1 : fib (2 * n + 2) = fib (2 * n) + fib (2 * n + 1) := rfl
    have e2 : fib (2 * n + 3) = fib (2 * n + 1) + fib (2 * n + 2) := rfl
    have k1 : 2 * (n + 1) = 2 * n + 2 := by omega
    rw [k1]
    show 2 * (ab n).1 + (ab n).2 = fib (2 * n + 2) ∧
      ((ab n).1 + (ab n).2) + (2 * (ab n).1 + (ab n).2) = fib (2 * n + 3)
    omega

/-- the automaton count agrees with the enumeration of R -/
theorem ab_counts_R :
    (List.range 7).all (fun i => ((stable (i + 1)).filter inR).length ==
      (ab (i + 1)).1 + (ab (i + 1)).2) = true := by decide +kernel

/-- relaxation time of Theorem C -/
def weight (s : List Nat) : Nat :=
  (List.range s.length).foldl (fun acc i => acc + (i + 1) * s.getD i 0) 0

def K (s : List Nat) : Nat :=
  let L := s.length
  L * (L + 1) - weight s + (if s.getLast? == some 2 then 1 else 0)

/-- Example ex:L2 -/
theorem K_L2 :
    [[0,0],[1,0],[2,0],[0,1],[1,1],[2,1],[0,2],[1,2],[2,2]].map K = [6,5,4,4,3,2,3,2,1] := by
  decide

/-- Lemma recgrains and Corollary recK, checked for L ≤ 7: on R the weight is at least
    L(L+1)/2 + [σ_L = 2], so K ≤ L(L+1)/2, with equality at (1,...,1). -/
def recKOK (L : Nat) : Bool :=
  ((stable L).filter inR).all (fun s =>
     weight s ≥ L * (L + 1) / 2 + (if s.getLast? == some 2 then 1 else 0) &&
     K s ≤ L * (L + 1) / 2) &&
  K (List.replicate L 1) == L * (L + 1) / 2

theorem recK_upto_7 : (List.range 7).all (fun i => recKOK (i + 1)) = true := by decide +kernel

/-- Corollary D(1): K ≤ L(L+1) with equality only at the empty pile, for L ≤ 7 -/
def maxKOK (L : Nat) : Bool :=
  (stable L).all (fun s => K s ≤ L * (L + 1) && (K s == L * (L + 1) → s.all (· == 0)))

theorem maxK_upto_7 : (List.range 7).all (fun i => maxKOK (i + 1)) = true := by decide +kernel
