"""
taxon_colors.py

Single source of truth for taxon colors across the pipeline, so every figure
(composition stacked bars in 06b, per-member effect plots in 08, and anything
added later) shows the same strain in the same color.

Palette: NatParksPalettes::natparks.pals(name = "DeathValley", n = 11,
type = "continuous") -- hex values below were computed with the package's own
colorRampPalette call, so they match R exactly.

Only the 11 real community members get palette colors. Unassigned reads and
the ambiguous Leaf59+p25C2 call are NOT community members, so they get
neutral grays rather than palette colors -- giving them palette colors would
visually present noise/ambiguity as if it were a strain. (This is also why
the palette is n = 11, not 13.)

Usage (scripts in the same directory):
    from taxon_colors import TAXON_COLORS, LOW_CONFIDENCE_COLOR
"""

# The 11 inoculated members, in the fixed order colors are assigned. Same
# order as 06b's TAXON_ORDER, so composition plots step through the palette
# in their stacking order.
MEMBER_ORDER = ['Leaf13', 'Leaf132', 'Leaf179', 'Leaf201', 'Leaf21', 'Leaf233',
                'Leaf234', 'Leaf391', 'Leaf59', 'Leaf84', 'p25C2']

# natparks.pals("DeathValley", n = 11, type = "continuous")
DEATHVALLEY_11 = ['#8C2B0E', '#AE5020', '#D07735', '#F2A450', '#9F7E59', '#132F5B',
                  '#2F4B7A', '#4A5982', '#60485B', '#865A63', '#B47E83']

TAXON_COLORS = dict(zip(MEMBER_ORDER, DEATHVALLEY_11))

# Non-member read categories (06b composition plots) -- neutral grays.
TAXON_COLORS['Unassigned'] = '#B0B0B0'
TAXON_COLORS['Leaf59+p25C2 (ambiguous)'] = '#5A5A5A'

# Fallback for any label not listed above (so a new label is visible, not a crash).
FALLBACK_COLOR = '#7F7F7F'

# Low-confidence taxa (low_confidence_flag in 07) are drawn in this color, at
# reduced opacity, INSTEAD of their taxon color in the effect plots -- a
# neutral gray that doesn't occur in the palette, so a flagged point can never
# be mistaken for a confident one of a similar hue.
LOW_CONFIDENCE_COLOR = '#9A9A9A'
LOW_CONFIDENCE_ALPHA = 0.45


def taxon_color(taxon):
    """Palette color for a taxon label, or the fallback gray."""
    return TAXON_COLORS.get(taxon, FALLBACK_COLOR)
