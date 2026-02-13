# UOR Foundation Research

Top-level container for UOR Foundation research projects. Each project lives as an independent repository and is included here as a [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) so the umbrella repo stays in sync with upstream while preserving each project’s own history and releases.

## Projects

| Project | Description |
|--------|-------------|
| [atlas-embeddings](./atlas-embeddings) | First-principles construction of exceptional Lie groups from the Atlas of Resonance Classes |

## Cloning this repository

To clone **research** and get all submodules (e.g. `atlas-embeddings`) in one step:

```bash
git clone --recurse-submodules https://github.com/UOR-Foundation/research.git
cd research
```

If you already cloned without submodules:

```bash
git submodule update --init --recursive
```

## Updating submodules

After pulling the latest **research** (in case the submodule pointers changed):

```bash
git pull
git submodule update --init --recursive
```

To update a specific project to its latest upstream commit and record that in **research**:

```bash
cd atlas-embeddings
git fetch origin
git checkout main
git pull origin main
cd ..
git add atlas-embeddings
git commit -m "Update atlas-embeddings submodule"
git push origin main
```

See [CONTRIBUTING.md](./CONTRIBUTING.md) for day-to-day submodule workflow and contributing to individual projects.

## License

MIT License. See [LICENSE](./LICENSE). Individual projects retain their own licenses.
