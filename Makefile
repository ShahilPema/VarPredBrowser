# VarPredBrowser Development Makefile

.PHONY: setup install install-python install-npm dev backend frontend build clean help

VENV := .venv
PYTHON := $(VENV)/bin/python
PIP := $(VENV)/bin/pip
FRONTEND_DIR := browser/frontend

# Default target
help:
	@echo "VarPredBrowser Development Commands"
	@echo ""
	@echo "  make setup          Full first-time setup (venv + all deps)"
	@echo "  make install        Install all dependencies"
	@echo "  make install-python Install Python dependencies only"
	@echo "  make install-npm    Install npm dependencies only"
	@echo "  make dev            Run backend + frontend dev servers"
	@echo "  make backend        Run backend server only"
	@echo "  make frontend       Run frontend dev server only"
	@echo "  make build          Build frontend for production"
	@echo "  make clean          Remove venv, node_modules, build artifacts"

# Create virtual environment
$(VENV)/bin/activate:
	python3 -m venv $(VENV)
	$(PIP) install --upgrade pip

# Full setup
setup: $(VENV)/bin/activate install
	@echo ""
	@echo "Setup complete! Run 'make dev' to start development servers."
	@echo "Activate venv with: source $(VENV)/bin/activate"

# Install all dependencies
install: install-python install-npm

# Install Python dependencies
install-python: $(VENV)/bin/activate
	$(PIP) install -r requirements-browser.txt

# Install npm dependencies
install-npm:
	cd $(FRONTEND_DIR) && npm install

# Run both backend and frontend (backend in background)
dev:
	@echo "Starting backend on http://localhost:8000"
	@echo "Starting frontend on http://localhost:5173"
	@echo "Press Ctrl+C to stop both servers"
	@trap 'kill 0' EXIT; \
	$(PYTHON) -m browser.backend.app & \
	cd $(FRONTEND_DIR) && npm run dev

# Run backend only
backend:
	$(PYTHON) -m browser.backend.app

# Run frontend dev server only
frontend:
	cd $(FRONTEND_DIR) && npm run dev

# Build frontend for production
build:
	cd $(FRONTEND_DIR) && npm run build

# Clean build artifacts and dependencies
clean:
	rm -rf $(VENV)
	rm -rf $(FRONTEND_DIR)/node_modules
	rm -rf $(FRONTEND_DIR)/dist
	@echo "Cleaned venv, node_modules, and dist"
