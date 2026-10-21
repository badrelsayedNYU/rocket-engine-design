# Kickoff Prompt (paste this into ChatGPT along with files 00–03)

---

I'm building a liquid rocket engine design tool in Python. I've attached four spec files:

- `00_PROJECT_BRIEF.md`: goals, principles, tech stack, repo layout, API
- `01_THEORY.md`: every equation to implement, with assumptions and references
- `02_STAGES.md`: the build plan, one stage at a time, with tests
- `03_VALIDATION.md`: real engine data to validate against (verify every number)

Please:
1. Read all four files fully before writing code.
2. Restate the plan in 5–8 bullets so I know you understood it, and flag anything in the spec that seems wrong or ambiguous.
3. Then build **Stage 0 only**, and stop.

Rules for the whole project:
- One stage at a time. Stop after each and wait for my go-ahead.
- Show complete files, not fragments, so I can paste them directly.
- Every physics function gets a docstring with the equation, its assumptions, and a reference.
- SI units internally. No unexplained constants.
- Verify any real-world number (engine data, CEA propellant names, Rao angle tables) instead of guessing, and tell me the source.
- After each stage, give me 2–3 sentences on the physics I should be able to explain in an interview.

---

## Tips for using this with ChatGPT

- If it starts writing several stages at once, tell it: *"Stop. Only Stage N. Follow 02_STAGES.md."*
- If a test fails, paste the full error back rather than describing it.
- Run each stage locally before moving on. Don't trust output you haven't run.
- If `rocketcea` fails to install (it needs a Fortran compiler on some systems), ask it to switch the `thermo.py` backend to Cantera; the brief allows this.
- Keep the files in the repo under `docs/spec/` so the AI (and anyone reading your repo) can see the design intent.
