// ResultScreen.tsx — shows prediction result

import React from "react";
import { View, Text, ScrollView, TouchableOpacity } from "react-native";
import { PredictionResult } from "../services/api";

export default function ResultScreen({ route, navigation }: any) {
  const { result, conditions }: {
    result: PredictionResult;
    conditions: any;
  } = route.params;

  const color = result.success ? "#1B5E20" : "#B71C1C";
  const bg    = result.success ? "#E8F5E9" : "#FFEBEE";

  return (
    <ScrollView style={{ flex: 1, backgroundColor: "#fff" }}
                contentContainerStyle={{ padding: 24 }}>

      {/* Reaction label */}
      <Text style={{ fontSize: 13, color: "#888", marginBottom: 4 }}>Reaction</Text>
      <Text style={{ fontWeight: "700", fontSize: 16, marginBottom: 20, color: "#333" }}>
        {conditions.amineType} ({conditions.position}-{conditions.substituent}) + DES
      </Text>

      {/* Success / Fail badge */}
      <View style={{
        backgroundColor: bg, borderRadius: 14,
        padding: 20, alignItems: "center", marginBottom: 24,
      }}>
        <Text style={{ fontSize: 40 }}>{result.success ? "✅" : "❌"}</Text>
        <Text style={{ fontSize: 22, fontWeight: "800", color, marginTop: 8 }}>
          {result.success ? "Reaction Likely" : "Reaction Unlikely"}
        </Text>
      </View>

      {/* Yield bar */}
      <Text style={{ fontWeight: "700", marginBottom: 6, color: "#333" }}>
        Predicted Yield: {result.yield_percent}%
      </Text>
      <View style={{ height: 14, backgroundColor: "#E0E0E0", borderRadius: 7, marginBottom: 20 }}>
        <View style={{
          height: 14,
          width: `${Math.min(result.yield_percent, 100)}%`,
          backgroundColor: color,
          borderRadius: 7,
        }} />
      </View>

      {/* Confidence */}
      <View style={{
        backgroundColor: "#F3F4F6", borderRadius: 10,
        padding: 12, marginBottom: 16,
      }}>
        <Text style={{ color: "#555" }}>
          <Text style={{ fontWeight: "700" }}>Confidence: </Text>
          {result.confidence}
        </Text>
      </View>

      {/* Warning */}
      {result.warning && (
        <View style={{
          backgroundColor: "#FFF8E1", borderRadius: 10,
          padding: 12, marginBottom: 16,
        }}>
          <Text style={{ color: "#F57F17" }}>{result.warning}</Text>
        </View>
      )}

      {/* Recommendation */}
      <View style={{
        backgroundColor: "#E3F2FD", borderRadius: 10,
        padding: 12, marginBottom: 24,
      }}>
        <Text style={{ fontWeight: "700", color: "#1565C0", marginBottom: 4 }}>
          Recommendation
        </Text>
        <Text style={{ color: "#1565C0" }}>{result.recommendation}</Text>
      </View>

      {/* Similar reactions from dataset */}
      {result.similar_reactions?.length > 0 && (
        <View style={{ marginBottom: 24 }}>
          <Text style={{ fontWeight: "700", color: "#333", marginBottom: 8 }}>
            Similar Reactions (from dataset)
          </Text>
          {result.similar_reactions.map((r, i) => (
            <View key={i} style={{
              flexDirection: "row", justifyContent: "space-between",
              paddingVertical: 6, borderBottomWidth: 1, borderColor: "#eee",
            }}>
              <Text style={{ color: "#555" }}>
                {r.amine} ({r.position}-{r.substituent})
              </Text>
              <Text style={{ fontWeight: "700", color: "#1565C0" }}>{r.yield}%</Text>
            </View>
          ))}
        </View>
      )}

      {/* Back button */}
      <TouchableOpacity
        onPress={() => navigation.goBack()}
        style={{
          backgroundColor: "#1565C0", padding: 16,
          borderRadius: 12, alignItems: "center",
        }}
      >
        <Text style={{ color: "#fff", fontWeight: "800" }}>Try Another Reaction</Text>
      </TouchableOpacity>

    </ScrollView>
  );
}
