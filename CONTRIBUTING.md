# Contributing to UOR Foundation Research

This repo is an **umbrella**: it only tracks which commit of each project to use. Development and contributions happen in the individual project repositories (e.g. [atlas-embeddings](https://github.com/UOR-Foundation/atlas-embeddings)).

## Working with submodules

### Clone the umbrella repo with all projects

```bash
git clone --recurse-submodules https://github.com/UOR-Foundation/research.git
cd research
```

You’ll have each project at the commit recorded by **research**.

### Work on a project (e.g. atlas-embeddings)

1. Go into the project directory and use it like a normal repo (branch, commit, push to **that** project’s GitHub repo):

   ```bash
   cd atlas-embeddings
   git checkout main
   git pull origin main
   # make changes, then:
   git add .
   git commit -m "Your change"
   git push origin main
   ```

2. From the **research** repo root, update the submodule pointer to the new commit and push:

   ```bash
   cd ..   # back to research/
   git add atlas-embeddings
   git commit -m "Update atlas-embeddings submodule"
   git push origin main
   ```

So: **code changes and PRs go to the project repo**; **research** only gets a commit that points to the new submodule commit.

### After someone else updates a submodule pointer

When you `git pull` in **research** and the submodule pointer changed:

```bash
git pull
git submodule update --init --recursive
```

That checks out the project at the commit now recorded by **research**.

### Pulling latest project code without changing the umbrella

To try the latest **atlas-embeddings** in your working tree without updating what **research** points to:

```bash
cd atlas-embeddings
git fetch origin
git checkout main
git pull origin main
```

Your **research** clone still points to the previous commit until you run `git add atlas-embeddings` and commit in **research**.

## Adding a new project as a submodule

From the **research** repo root:

```bash
git submodule add https://github.com/UOR-Foundation/<project>.git <project>
git add .gitmodules <project>
git commit -m "Add <project> as submodule"
git push origin main
```

## Quick reference

| Task | Where | Command |
|------|--------|--------|
| Clone research with all projects | Anywhere | `git clone --recurse-submodules https://github.com/UOR-Foundation/research.git` |
| Init submodules after clone | Inside `research/` | `git submodule update --init --recursive` |
| Sync to recorded commits after pull | Inside `research/` | `git submodule update --init --recursive` |
| Develop / contribute | Inside `research/atlas-embeddings/` | Normal git: branch, commit, push to **atlas-embeddings** repo |
| Record new project commit in research | Inside `research/` | `git add atlas-embeddings` then commit and push |

Contributions (issues, PRs, docs) for each project should go to that project’s repository, not to this umbrella repo.
