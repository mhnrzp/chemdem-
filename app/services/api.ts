// api.ts — all HTTP calls to the Python backend

const BASE_URL = "http://localhost:8000"; // change to your machine IP when testing on real device

export interface ReactionRequest {
  amine_type:  string;
  substituent: string;
  position:    string;
  temperature: string;
  catalyst:    boolean;
}

export interface SimilarReaction {
  amine: string;
  substituent: string;
  position: string;
  yield: number;
}

export interface PredictionResult {
  success:           boolean;
  yield_percent:     number;
  confidence:        string;
  message:           string;
  warning:           string | null;
  recommendation:    string;
  similar_reactions: SimilarReaction[];
}

export async function predictReaction(req: ReactionRequest): Promise<PredictionResult> {
  const res = await fetch(`${BASE_URL}/predict`, {
    method:  "POST",
    headers: { "Content-Type": "application/json" },
    body:    JSON.stringify(req),
  });
  if (!res.ok) throw new Error(`Server error: ${res.status}`);
  return res.json();
}
