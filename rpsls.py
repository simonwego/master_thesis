"""
Simulerer evolusjon av en populasjon som spiller RPSLS (rock-paper-scissors-lizard-spock)
i dyadiske interaksjoner, med gruppevise (typevise) og justerbare biaser for valg av trekk.

Idé:
- Populasjonen består av "grupper" (typer/strategier).
- Hver gruppe har en basis-bias (logitter) for de 5 trekkene.
- I tillegg kan du spesifisere dyadiske bias-justeringer: når gruppe g møter h,
  får g et ekstra logit-tillegg (vektor) på sine trekk-sannsynligheter.
- I hver generasjon:
  1) Simuler mange kamper (Monte Carlo) mellom individer trukket fra populasjonen.
  2) Beregn gruppevis gjennomsnitts-payoff (fitness).
  3) Oppdater gruppefrekvenser med en Wright–Fisher-lignende seleksjon.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import numpy as np

# ---------- RPSLS payoff ----------

MOVES = ["rock", "paper", "scissors", "lizard", "spock"]
M = len(MOVES)

# payoff[a, b] = payoff til spiller A som spiller move a mot B som spiller move b
# +1 hvis A vinner, -1 hvis A taper, 0 hvis uavgjort
def rpsls_payoff_matrix() -> np.ndarray:
    idx = {m: i for i, m in enumerate(MOVES)}
    wins_against = {
        "rock":     ["scissors", "lizard"],
        "paper":    ["rock", "spock"],
        "scissors": ["paper", "lizard"],
        "lizard":   ["paper", "spock"],
        "spock":    ["rock", "scissors"],
    }
    P = np.zeros((M, M), dtype=np.int8)
    for a in MOVES:
        for b in MOVES:
            ia, ib = idx[a], idx[b]
            if a == b:
                P[ia, ib] = 0
            elif b in wins_against[a]:
                P[ia, ib] = 1
            else:
                P[ia, ib] = -1
    return P


PAYOFF = rpsls_payoff_matrix()


# ---------- Utility ----------

def softmax(logits: np.ndarray) -> np.ndarray:
    x = np.asarray(logits, dtype=float)
    x = x - np.max(x)
    ex = np.exp(x)
    return ex / np.sum(ex)

def clamp_min(x: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    return np.maximum(x, eps)


# ---------- Model ----------

@dataclass(frozen=True)
class GroupSpec:
    name: str
    base_logits: np.ndarray  # shape (5,)

@dataclass
class SimulatorConfig:
    selection_strength: float = 1.0   # høyere => sterkere seleksjon
    matches_per_generation: int = 50_000
    # hvis True: dyadisk bias brukes, hvis False: kun base_logits
    use_pair_bias: bool = True
    # små "mutasjoner"/drift i frekvenser kan legges inn som en Dirichlet-smoothing
    dirichlet_smoothing: float = 0.0  # f.eks 1e-3 for litt smoothing

class RPSLSEvolutionSimulator:
    def __init__(
        self,
        group_specs: List[GroupSpec],
        initial_counts: np.ndarray,
        pair_bias_logits: Optional[Dict[Tuple[int, int], np.ndarray]] = None,
        seed: Optional[int] = None,
    ):
        """
        pair_bias_logits[(g, h)] = logit-tillegg (shape (5,)) som gruppe g får når den møter gruppe h.
        Hvis en (g,h) ikke finnes: 0-tillegg.
        """
        self.rng = np.random.default_rng(seed)

        self.groups = group_specs
        self.G = len(group_specs)
        self.counts = np.array(initial_counts, dtype=int)
        if self.counts.shape != (self.G,):
            raise ValueError(f"initial_counts må ha shape ({self.G},), fikk {self.counts.shape}")

        self.base_logits = np.stack([g.base_logits for g in group_specs], axis=0)  # (G, 5)
        if self.base_logits.shape != (self.G, M):
            raise ValueError(f"base_logits må ha shape (G,5). Fikk {self.base_logits.shape}")

        self.pair_bias = pair_bias_logits or {}

    def total_pop(self) -> int:
        return int(self.counts.sum())

    def group_probs(self) -> np.ndarray:
        n = self.total_pop()
        if n <= 0:
            raise ValueError("Populasjonsstørrelse må være > 0")
        return self.counts / n

    def _logits_for(self, g: int, h: int, use_pair_bias: bool) -> np.ndarray:
        logits = self.base_logits[g].copy()
        if use_pair_bias:
            add = self.pair_bias.get((g, h), None)
            if add is not None:
                logits = logits + add
        return logits

    def _sample_group_pairings(self, n_pairs: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Trekker n_pairs dyadiske møter (gruppeindekser) fra populasjonsfordelingen.
        For enkelhet: trekker uavhengig med tilbakelegging på gruppenivå.
        """
        p = self.group_probs()
        gA = self.rng.choice(self.G, size=n_pairs, p=p)
        gB = self.rng.choice(self.G, size=n_pairs, p=p)
        return gA, gB

    def _simulate_generation_payoffs(self, cfg: SimulatorConfig) -> np.ndarray:
        """
        Returnerer estimert gjennomsnittlig payoff (fitness) per gruppe, shape (G,).
        """
        n_pairs = cfg.matches_per_generation
        gA, gB = self._sample_group_pairings(n_pairs)

        # For å være rask: batch per unik (g,h)
        payoff_sum = np.zeros(self.G, dtype=float)
        games_count = np.zeros(self.G, dtype=float)

        # Finn alle unike matchups
        pairs = np.stack([gA, gB], axis=1)
        uniq, inv = np.unique(pairs, axis=0, return_inverse=True)

        for k, (g, h) in enumerate(uniq):
            mask = (inv == k)
            m = int(mask.sum())
            if m == 0:
                continue

            p_g = softmax(self._logits_for(int(g), int(h), cfg.use_pair_bias))
            p_h = softmax(self._logits_for(int(h), int(g), cfg.use_pair_bias))

            # Sample trekk
            a_moves = self.rng.choice(M, size=m, p=p_g)
            b_moves = self.rng.choice(M, size=m, p=p_h)

            # Payoff til A (gruppe g) og til B (gruppe h)
            a_pay = PAYOFF[a_moves, b_moves].astype(float)
            b_pay = -a_pay  # zero-sum

            payoff_sum[int(g)] += a_pay.sum()
            payoff_sum[int(h)] += b_pay.sum()
            games_count[int(g)] += m
            games_count[int(h)] += m

        # Gjennomsnittlig payoff per match for hver gruppe
        fitness = payoff_sum / clamp_min(games_count)
        return fitness

    def step(self, cfg: SimulatorConfig) -> Dict[str, np.ndarray]:
        """
        Kjør én generasjon: simuler kamper -> beregn fitness -> oppdater counts.
        Returnerer diagnostikk.
        """
        N = self.total_pop()
        fitness = self._simulate_generation_payoffs(cfg)

        # Konverter payoff til reproduktiv vekt:
        # w_g ∝ exp(s * fitness_g)
        s = float(cfg.selection_strength)
        w = np.exp(s * fitness)
        w = w / w.sum()

        # Wright–Fisher oppdatering (multinomial)
        new_counts = self.rng.multinomial(N, w)

        # Valgfri smoothing (dirichlet-ish) for å unngå at typer dør helt pga sampling-støy
        if cfg.dirichlet_smoothing > 0:
            alpha = cfg.dirichlet_smoothing
            # bland litt med uniform
            u = np.full(self.G, 1.0 / self.G)
            w2 = (1 - alpha) * (new_counts / N) + alpha * u
            w2 = w2 / w2.sum()
            new_counts = self.rng.multinomial(N, w2)

        self.counts = new_counts

        return {
            "fitness": fitness,
            "repro_probs": w,
            "counts": new_counts.copy(),
            "freqs": new_counts / N,
        }

    def run(self, generations: int, cfg: SimulatorConfig) -> Dict[str, np.ndarray]:
        """
        Kjør flere generasjoner og lagrer historikk.
        """
        G = self.G
        counts_hist = np.zeros((generations + 1, G), dtype=int)
        freqs_hist = np.zeros((generations + 1, G), dtype=float)
        fitness_hist = np.zeros((generations, G), dtype=float)

        counts_hist[0] = self.counts
        freqs_hist[0] = self.group_probs()

        for t in range(generations):
            out = self.step(cfg)
            fitness_hist[t] = out["fitness"]
            counts_hist[t + 1] = out["counts"]
            freqs_hist[t + 1] = out["freqs"]

        return {
            "counts": counts_hist,
            "freqs": freqs_hist,
            "fitness": fitness_hist,
            "group_names": np.array([g.name for g in self.groups], dtype=object),
        }


# ---------- Eksempelbruk ----------

if __name__ == "__main__":
    # Tre grupper med forskjellige "preferanser" (base_logits)
    group_specs = [
        GroupSpec("Rocky",    np.array([2.0, 0.0, 0.0, 0.0, 0.0])),  # liker rock
        GroupSpec("Papery",   np.array([0.0, 2.0, 0.0, 0.0, 0.0])),  # liker paper
        GroupSpec("Spocky",   np.array([0.0, 0.0, 0.0, 0.0, 2.0])),  # liker spock
    ]

    N = 3000
    initial_counts = np.array([1000, 1000, 1000])

    # Dyadiske biaser: gruppe g får ekstra logit-tillegg når den møter h.
    # Her: "Rocky" (0) blir ekstra glad i lizard når den møter "Spocky" (2),
    # og "Papery" (1) blir ekstra glad i scissors når den møter "Rocky" (0).
    pair_bias = {
        (0, 2): np.array([0.0, 0.0, 0.0, 1.5, 0.0]),  # vs Spocky => mer lizard
        (1, 0): np.array([0.0, 0.0, 1.5, 0.0, 0.0]),  # vs Rocky  => mer scissors
    }

    sim = RPSLSEvolutionSimulator(
        group_specs=group_specs,
        initial_counts=initial_counts,
        pair_bias_logits=pair_bias,
        seed=123,
    )

    cfg = SimulatorConfig(
        selection_strength=2.0,
        matches_per_generation=80_000,
        use_pair_bias=True,
        dirichlet_smoothing=0.0,
    )

    hist = sim.run(generations=60, cfg=cfg)

    # Print siste frekvenser
    names = hist["group_names"]
    final_freqs = hist["freqs"][-1]
    print("Final frequencies:")
    for n, f in zip(names, final_freqs):
        print(f"  {n:>8s}: {f:.3f}")

    # (Valgfritt) Plot hvis du vil:
    try:
        import matplotlib.pyplot as plt
        t = np.arange(hist["freqs"].shape[0])
        for i, n in enumerate(names):
            plt.plot(t, hist["freqs"][:, i], label=str(n))
        plt.xlabel("Generation")
        plt.ylabel("Frequency")
        plt.legend()
        plt.tight_layout()
        plt.show()
    except Exception:
        pass