import numpy as np

# --- RPSLS-regler ---
# Vi koder trekk som: 0=rock, 1=paper, 2=scissors, 3=lizard, 4=spock
# "beats[a]" er settet av trekk a slår.
beats = {
    0: {2, 3},  # rock beats scissors, lizard
    1: {0, 4},  # paper beats rock, spock
    2: {1, 3},  # scissors beats paper, lizard
    3: {1, 4},  # lizard beats paper, spock
    4: {0, 2},  # spock beats rock, scissors
}

def round_outcome(a: int, b: int) -> int:
    """
    Returnerer  1 hvis a vinner, -1 hvis b vinner, 0 hvis uavgjort.
    """
    if a == b:
        return 0
    if b in beats[a]:
        return 1
    return -1

def play_series(rng: np.random.Generator, wins_to_take: int = 5) -> int:
    """
    Simulerer en serie der to spillere velger uniformt.
    Returnerer 0 hvis spiller A vinner serien, 1 hvis spiller B vinner serien.
    Uavgjorte runder teller ikke (spilles om).
    """
    wA = 0
    wB = 0
    while wA < wins_to_take and wB < wins_to_take:
        a = int(rng.integers(0, 5))
        b = int(rng.integers(0, 5))
        o = round_outcome(a, b)
        if o == 1:
            wA += 1
        elif o == -1:
            wB += 1
    return 0 if wA == wins_to_take else 1

def simulate_tournament(num_players: int, series_per_pair: int = 1000, wins_to_take: int = 5, seed: int = 0):
    """
    Returnerer:
      wins[i,j] = antall serier i spiller i vinner mot spiller j
      (wins[j,i] = series_per_pair - wins[i,j])
    """
    rng = np.random.default_rng(seed)
    wins = np.zeros((num_players, num_players), dtype=int)

    for i in range(num_players):
        for j in range(i + 1, num_players):
            w_i = 0
            for _ in range(series_per_pair):
                winner = play_series(rng, wins_to_take=wins_to_take)
                if winner == 0:
                    w_i += 1
            wins[i, j] = w_i
            wins[j, i] = series_per_pair - w_i

    return wins

# --- Eksempelbruk ---
if __name__ == "__main__":
    n_players = 8
    wins = simulate_tournament(num_players=n_players, series_per_pair=2000, wins_to_take=5, seed=42)

    print("wins[i,j] = antall serier i slo j")
    print(wins)
    print("\nSeiersandel (rad i mot kolonne j):")
    print(np.round(wins / (wins + wins.T + np.eye(n_players)), 3))  # litt penere, diagonal blir NaN-ish håndtert av +I