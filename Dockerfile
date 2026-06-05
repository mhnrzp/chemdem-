FROM python:3.11-slim

# Install system dependencies needed by RDKit
RUN apt-get update && apt-get install -y \
    libxrender1 \
    libxext6 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy and install Python dependencies
COPY backend/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy backend source
COPY backend/ .

# Expose port (7860 for HF Spaces, 8000 for Railway via PORT env var)
EXPOSE 7860

# Start FastAPI — use $PORT if set (Railway sets 8000), otherwise 7860 (HF Spaces default)
CMD ["sh", "-c", "uvicorn main:app --host 0.0.0.0 --port ${PORT:-7860}"]
