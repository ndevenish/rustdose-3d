# RADDOSE-3D Browser Demo

Interactive dose modelling in the browser via WebAssembly.

## Prerequisites

Install the Rust toolchain and the following tools if you don't have them:

```sh
# Rust (stable)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# wasm32 compile target
rustup target add wasm32-unknown-unknown

# wasm-pack (builds and packages the WASM module)
cargo install wasm-pack

# wasm-opt (Binaryen optimiser — used by wasm-pack release builds)
cargo install wasm-opt
```

## Build

Run from the `crates/raddose3d-wasm/` directory:

```sh
cd crates/raddose3d-wasm
wasm-pack build --target web --out-dir pkg
```

This produces `pkg/` containing:

| File | Description |
|------|-------------|
| `raddose3d_wasm_bg.wasm` | Compiled + optimised WebAssembly binary |
| `raddose3d_wasm.js` | ES module loader / JS glue |
| `raddose3d_wasm.d.ts` | TypeScript type declarations |
| `package.json` | npm package metadata |

The `web/index.html` demo imports from `../pkg/` (relative to itself), so the
directory structure must be:

```
crates/raddose3d-wasm/
├── pkg/
│   ├── raddose3d_wasm_bg.wasm
│   ├── raddose3d_wasm.js
│   └── ...
└── web/
    └── index.html       ← serve this
```

## Serve

The page **must be served over HTTP** — browsers block ES module imports and
WebAssembly instantiation from `file://` URLs.

Any static file server works. The server root should be `crates/raddose3d-wasm/`
so that both `web/` and `pkg/` are reachable.

### Python (no install required)

```sh
cd crates/raddose3d-wasm
python3 -m http.server 8080
```

Then open <http://localhost:8080/web/>.

### Node.js (`npx serve`)

```sh
cd crates/raddose3d-wasm
npx serve .
```

### Other options

```sh
# cargo-server (Rust)
cargo install cargo-server
cargo server --open

# miniserve (Rust)
cargo install miniserve
miniserve . --index web/index.html
```

## Tests

Tests run the full simulation inside Node.js via `wasm-bindgen-test`:

```sh
cd crates/raddose3d-wasm
wasm-pack test --node
```
