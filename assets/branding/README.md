# rsx branding assets

Primary logos and icons for the rsx (rsx-rs) project.

## Files

### Recommended (vector, preferred)

- `rsx-logo.svg` — full horizontal logotype (icon + "rsx" wordmark). Use in READMEs, docs headers, website.
- `rsx-icon.svg` — square icon mark only. Use for favicons, GitHub avatar, social icons, app icons, PyPI.
- `rsx-logo-dark.svg` / `rsx-icon-dark.svg` — light-on-dark variants for dark mode headers, GitHub dark theme, etc.

### PNG fallbacks (transparent background)

- `rsx-logo.png` / `rsx-logo@2x.png`
- `rsx-icon.png` / `rsx-icon@2x.png`
- `...-dark...` equivalents

### Concepts

The `concepts/` directory contains earlier generative explorations (marker streams, DNA helix interpretations, monograms). Use only for reference or inspiration; the SVGs above are the canonical assets.

## Colors (official)

- Teal (primary): `#0F766E` (light) / `#14B8A6` (dark variant)
- Orange accent (significant markers / highlights): `#EA580C` (light) / `#F97316` (dark)
- Slate (text/cut): `#0F172A` (light) / `#F1F5F9` (dark)

## Usage

### README header (typical)

```markdown
<p align="center">
  <img src="assets/branding/rsx-logo.svg" alt="rsx" width="420">
</p>
```

Or for GitHub that sometimes prefers PNG:

```markdown
<p align="center">
  <img src="assets/branding/rsx-logo.png" alt="rsx" width="420">
</p>
```

### Sphinx docs (shibuya or alabaster)

In `docs/source/conf.py`:

```python
html_logo = "../assets/branding/rsx-icon.svg"
# or html_favicon = "../assets/branding/rsx-icon.svg"
```

### Favicon

Use a 32x32 or 16x16 export from `rsx-icon.png`, or run through https://realfavicongenerator.net on the SVG/PNG.

### Dark mode

Serve `rsx-logo-dark.svg` when the user's `prefers-color-scheme: dark` (or use CSS to invert/recolor the light asset).

## License

These assets are part of the rsx-rs project and are licensed under the same terms as the software (GPL-3.0-or-later). Attribution appreciated but not required for normal use.

