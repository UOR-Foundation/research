# Contributing to UOR Foundation Research

This repo is a **monorepo**: each project is a subdirectory (e.g. `atlas-embeddings/`). Contribute by editing files in the relevant project and opening a pull request against **research**.

## Working on a project

1. Clone the repo and go into the project directory:

   ```bash
   git clone https://github.com/UOR-Foundation/research.git
   cd research/atlas-embeddings
   ```

2. Make changes, run tests, and commit from the repo root (so history stays in research):

   ```bash
   cd research
   git add atlas-embeddings/
   git commit -m "atlas-embeddings: your change description"
   git push origin main
   ```

3. Open a pull request on [UOR-Foundation/research](https://github.com/UOR-Foundation/research).

## Project layout

- **atlas-embeddings/** — Rust crate + Lean 4 formalization. See [atlas-embeddings/README.md](./atlas-embeddings/README.md) for build and test instructions.

Contributions (issues, PRs) go to this repository; specify the project (e.g. `atlas-embeddings`) in the PR title or description.
