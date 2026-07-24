import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import tailwindcss from '@tailwindcss/vite'

// For GitHub Pages PROJECT sites the app is served from a subpath
// (https://<user>.github.io/<repo>/), so the asset base must be that subpath.
// Override at build time:  VITE_BASE=/rmc-phonon-dynamics/ npm run build
// Defaults to '/' for local dev and user/organization Pages.
// https://vite.dev/config/
export default defineConfig({
  base: process.env.VITE_BASE || '/',
  plugins: [
    react(),
    tailwindcss(),
  ],
})
