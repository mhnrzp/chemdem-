// HomeScreen.tsx
// User picks reaction conditions → hits Predict → goes to ResultScreen

import React, { useState } from "react";
import {
  View, Text, TouchableOpacity, ScrollView, Switch
} from "react-native";
import { predictReaction } from "../services/api";

// ── Option lists (mirror backend /options) ───────────────────────────────────
const AMINE_TYPES   = ["aniline", "benzylamine"];
const SUBSTITUENTS  = ["none", "F", "Cl", "Br", "OH", "OMe", "CF3", "tolyl"];
const POSITIONS     = ["ortho", "meta", "para", "none"];
const TEMPERATURES  = ["r.t.", "reflux"];

// ── Small selector component ─────────────────────────────────────────────────
function Selector({ label, options, value, onChange }: {
  label: string;
  options: string[];
  value: string;
  onChange: (v: string) => void;
}) {
  return (
    <View style={{ marginBottom: 16 }}>
      <Text style={{ fontWeight: "700", marginBottom: 6, color: "#333" }}>{label}</Text>
      <View style={{ flexDirection: "row", flexWrap: "wrap", gap: 8 }}>
        {options.map((opt) => (
          <TouchableOpacity
            key={opt}
            onPress={() => onChange(opt)}
            style={{
              paddingVertical: 6,
              paddingHorizontal: 14,
              borderRadius: 20,
              backgroundColor: value === opt ? "#1565C0" : "#E3F2FD",
            }}
          >
            <Text style={{ color: value === opt ? "#fff" : "#1565C0", fontWeight: "600" }}>
              {opt}
            </Text>
          </TouchableOpacity>
        ))}
      </View>
    </View>
  );
}

// ── Main Screen ───────────────────────────────────────────────────────────────
export default function HomeScreen({ navigation }: any) {
  const [amineType,   setAmineType]   = useState("aniline");
  const [substituent, setSubstituent] = useState("none");
  const [position,    setPosition]    = useState("para");
  const [temperature, setTemperature] = useState("r.t.");
  const [catalyst,    setCatalyst]    = useState(false);
  const [loading,     setLoading]     = useState(false);
  const [error,       setError]       = useState<string | null>(null);

  async function handlePredict() {
    setLoading(true);
    setError(null);
    try {
      const result = await predictReaction({
        amine_type:  amineType,
        substituent,
        position,
        temperature,
        catalyst,
      });
      navigation.navigate("Result", { result, conditions: {
        amineType, substituent, position, temperature, catalyst
      }});
    } catch (e: any) {
      setError("Could not reach server. Is the backend running?");
    } finally {
      setLoading(false);
    }
  }

  return (
    <ScrollView style={{ flex: 1, backgroundColor: "#fff" }}
                contentContainerStyle={{ padding: 24 }}>

      {/* Header */}
      <Text style={{ fontSize: 26, fontWeight: "800", color: "#1565C0" }}>Chemdem</Text>
      <Text style={{ color: "#888", marginBottom: 24 }}>
        Squaric Acid Monoamide Synthesis Predictor
      </Text>

      {/* Reaction label */}
      <View style={{
        backgroundColor: "#F3F4F6", borderRadius: 10,
        padding: 12, marginBottom: 24
      }}>
        <Text style={{ fontFamily: "monospace", color: "#333", fontSize: 13 }}>
          {amineType} ({position !== "none" ? position + "-" : ""}{substituent}) +
          {"\n"}Diethyl Squarate (DES)
          {"\n"}→ Monosquarate-amide
        </Text>
      </View>

      {/* Selectors */}
      <Selector label="Amine Type"    options={AMINE_TYPES}  value={amineType}   onChange={setAmineType}   />
      <Selector label="Substituent"   options={SUBSTITUENTS} value={substituent} onChange={setSubstituent} />
      <Selector label="Position"      options={POSITIONS}    value={position}    onChange={setPosition}    />
      <Selector label="Temperature"   options={TEMPERATURES} value={temperature} onChange={setTemperature} />

      {/* Catalyst toggle */}
      <View style={{ flexDirection: "row", alignItems: "center", marginBottom: 28 }}>
        <Text style={{ fontWeight: "700", color: "#333", flex: 1 }}>
          Zn(OTf)₂ Catalyst
        </Text>
        <Switch value={catalyst} onValueChange={setCatalyst}
                trackColor={{ true: "#1565C0" }} />
        <Text style={{ marginLeft: 8, color: "#888" }}>
          {catalyst ? "Yes" : "No"}
        </Text>
      </View>

      {/* Error */}
      {error && (
        <Text style={{ color: "#C62828", marginBottom: 12 }}>{error}</Text>
      )}

      {/* Predict button */}
      <TouchableOpacity
        onPress={handlePredict}
        disabled={loading}
        style={{
          backgroundColor: loading ? "#90CAF9" : "#1565C0",
          padding: 16, borderRadius: 12, alignItems: "center",
        }}
      >
        <Text style={{ color: "#fff", fontWeight: "800", fontSize: 16 }}>
          {loading ? "Predicting..." : "Predict Outcome"}
        </Text>
      </TouchableOpacity>

    </ScrollView>
  );
}
