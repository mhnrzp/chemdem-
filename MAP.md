# Chemdem — Project Map

## What This App Does
Mobile app that predicts chemical reaction outcomes (success/fail + estimated yield).
University project. Uses IBM RXN AI API for predictions.

---

## Folder Structure

```
Chemdem/
│
├── backend/                    ← Python API (run this on your computer)
│   ├── main.py                 ← FastAPI app, defines /predict and /reactions endpoints
│   ├── predictor.py            ← IBM RXN API wrapper (swap here for other models)
│   └── requirements.txt        ← pip install -r requirements.txt
│
├── app/                        ← Expo (React Native) mobile app
│   ├── screens/
│   │   ├── HomeScreen.tsx      ← Input screen (reaction list + advanced SMILES toggle)
│   │   └── ResultScreen.tsx    ← Results screen (badge, yield bar, confidence)
│   ├── services/
│   │   └── api.ts              ← All HTTP calls to the backend
│   └── components/             ← (empty) shared UI components go here
│
└── design/
    └── figma-exports/          ← Drop Figma PNG/SVG exports here
```

---

## Data Flow

```
User picks reaction (or types SMILES)
        ↓
HomeScreen.tsx  →  POST /predict  →  main.py
                                        ↓
                                   predictor.py
                                        ↓
                                   IBM RXN API (free, online)
                                        ↓
                              { success, yield_percent, confidence }
        ↓
ResultScreen.tsx displays result
```

---

## API Endpoints

| Method | Path         | What it does                          |
|--------|--------------|---------------------------------------|
| POST   | /predict     | Takes SMILES + type, returns outcome  |
| GET    | /reactions   | Returns the built-in reaction list    |
| GET    | /health      | Sanity check                          |

---

## Built-in Reaction List

| ID             | Name                  | Category      |
|----------------|-----------------------|---------------|
| suzuki         | Suzuki Coupling       | coupling      |
| heck           | Heck Reaction         | coupling      |
| negishi        | Negishi Coupling      | coupling      |
| sn1            | SN1 Substitution      | substitution  |
| sn2            | SN2 Substitution      | substitution  |
| aldol          | Aldol Condensation    | condensation  |
| esterification | Esterification        | condensation  |
| grignard       | Grignard Reaction     | addition      |

---

## Setup (when ready to run)

### Backend
```bash
cd backend
pip install -r requirements.txt
# Add your IBM RXN API key to predictor.py
uvicorn main:app --reload
# Runs at http://localhost:8000
```

### Mobile App
```bash
npm install -g expo-cli
npx create-expo-app chemdem --template blank-typescript
# Copy screens/ services/ into the new Expo project
npx expo start
```

---

## Current Status

- [x] Folder structure mapped
- [x] Backend skeleton (main.py, predictor.py)
- [x] App screens skeleton (HomeScreen, ResultScreen)
- [x] API service layer (api.ts)
- [ ] IBM RXN API key added
- [ ] Figma design created
- [ ] Figma styles applied to screens
- [ ] Full integration test
