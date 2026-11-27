#!/usr/bin/env python3
"""
Purpose:
    Composite Figure 1 (`Figure1_panel_a_updated.png`, `Figure1_panel_b_cleaned.png`, combined
    Figure 1, and clean website UI crops) from high-resolution screenshots of
    https://deepcre.ipk-gatersleben.de/ and the cleaned Figure 1b schematic (`No. of genes`,
    `log10(mean TPM+1)`, `gDR`).

Manuscript (Peleke et al., 2026):
    Builds Figure 1 (panels a & b). Note: reads external screenshot and crop PNGs from
    `EXTERNAL/figure_inputs/` (`fig1_ab_crop.png`, `deepcre_random100_*.png`,
    `tutorial_canvas_10.png`) that are not stored in this repository.

Usage & I/O (run from repo root with Pillow + matplotlib installed and `EXTERNAL/figure_inputs/*.png` present):
    python3 code/04_figures/build_updated_fig1.py
    - Writes cropped panels and composites into `EXTERNAL/figure_output/figures/`

How to change or extend:
    Update `FIG1_INPUT_DIR`, `OUT_DIR`, font paths (`FONT_REG`, `FONT_BOLD`), or screenshot crop coordinates if
    re-rendering against a newer web UI capture.
"""
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageOps
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIG1_INPUT_DIR = Path("EXTERNAL/figure_inputs")  # EXTERNAL (figure_inputs): see EXTERNAL_INPUTS.md
OUT_DIR = Path("EXTERNAL/figure_output/figures")  # EXTERNAL (figure_output): see EXTERNAL_INPUTS.md
OUT_DIR.mkdir(parents=True, exist_ok=True)

FONT_REG = "/usr/share/fonts/opentype/urw-base35/URWGothic-Book.otf"
FONT_BOLD = "/usr/share/fonts/opentype/urw-base35/URWGothic-Demi.otf"


def make_clean_panel_b() -> Image.Image:
    im = Image.open(FIG1_INPUT_DIR / "fig1_ab_crop.png").convert("RGB")
    arr = np.array(im).copy()

    # 1. Remove leaf, chromosome, chromatin beads, and "No.# of genes" to the left of x=1332
    arr[120:402, 1040:1332] = 255
    arr[380:440, 1040:1250] = 255
    # 2. Remove floating RNA wave icon on the right (x=1680..1870, y=220..335)
    arr[220:335, 1680:1870] = 255
    # 3. Replace 'd' in dDR with 'g' from gUR
    arr[804:840, 1468:1488] = arr[804:840, 1084:1104]

    im_b = Image.fromarray(arr)
    draw = ImageDraw.Draw(im_b)

    # 4. Replace log2Max(TPM+1) with log10(mean TPM+1)
    draw.rectangle([1355, 415, 1630, 464], fill=(255, 255, 255))
    font_axis = ImageFont.truetype(FONT_REG, 33)
    tmp = Image.new("RGB", (340, 52), (255, 255, 255))
    tdraw = ImageDraw.Draw(tmp)
    tdraw.text((0, 6), "log10(mean TPM+1)", fill=(54, 53, 53), font=font_axis)
    tmp_scaled = tmp.resize((int(340 * 0.90), 52), resample=Image.BICUBIC)
    im_b.paste(tmp_scaled, (1358, 415))

    # 5. Draw clean "No. of\ngenes" right-aligned to the left of the histogram y-axis
    font_y = ImageFont.truetype(FONT_REG, 34)
    tmp_y = Image.new("RGB", (165, 110), (255, 255, 255))
    ydraw = ImageDraw.Draw(tmp_y)
    for idx, line in enumerate(["No. of", "genes"]):
        w, h = ydraw.textsize(line, font=font_y)
        ydraw.text((158 - w, 6 + idx * 44), line, fill=(54, 53, 53), font=font_y)
    tmp_y_scaled = tmp_y.resize((int(165 * 0.92), 110), resample=Image.BICUBIC)
    im_b.paste(tmp_y_scaled, (1325 - tmp_y_scaled.width, 255))

    # Also erase the old left curly brace (x=900..1040) so we can place a clean divider or brace
    arr_b = np.array(im_b)
    arr_b[:, 900:1035] = 255
    im_b_clean = Image.fromarray(arr_b).crop((1020, 0, 1925, 1470))
    im_b_clean.save(OUT_DIR / "Figure1_panel_b_cleaned.png", dpi=(300, 300))
    return im_b_clean


def add_card_border(im: Image.Image, border_color=(215, 212, 225), pad=6) -> Image.Image:
    im_pad = ImageOps.expand(im, border=pad, fill=(255, 255, 255))
    im_brd = ImageOps.expand(im_pad, border=2, fill=border_color)
    return im_brd


def make_updated_panel_a(target_height: int = 1470) -> Image.Image:
    """
    Create a crisp, high-resolution Panel a (width 1120 px, height 1470 px)
    that guides the user through the live website in 3 clear steps:
      1. Inputs (Sidebar screenshot + clear step callouts, consistent A. thaliana example)
      2. Processing & Modelling on web server
      3. Interactive Outputs (Predictions, Explanation, Mutation controls + result)
    """
    W, H = 1120, target_height
    canvas = Image.new("RGB", (W, H), (255, 255, 255))
    draw = ImageDraw.Draw(canvas)

    f_panel = ImageFont.truetype(FONT_BOLD, 50)
    f_title = ImageFont.truetype(FONT_BOLD, 30)
    f_sec = ImageFont.truetype(FONT_BOLD, 27)
    f_sub = ImageFont.truetype(FONT_BOLD, 24)
    f_body = ImageFont.truetype(FONT_REG, 21)
    f_small = ImageFont.truetype(FONT_REG, 17)
    f_step = ImageFont.truetype(FONT_BOLD, 20)

    PURPLE = (88, 44, 131)
    DARK = (35, 31, 32)
    GREY = (115, 115, 120)
    LIGHT_BG = (248, 246, 252)
    STEP_BG = (88, 44, 131)

    # Panel letter 'a'
    draw.text((18, 12), "a", fill=(0, 0, 0), font=f_panel)
    draw.text((22, 82), "deepCRE", fill=PURPLE, font=f_title)
    w_dc, _ = draw.textsize("deepCRE ", font=f_title)
    draw.text((22 + w_dc, 82), "web toolkit workflow", fill=DARK, font=f_title)

    # Crop and save individual live website screenshots
    im_pred_full = Image.open(FIG1_INPUT_DIR / "deepcre_random100_predictions.png").convert("RGB")
    im_expl_full = Image.open(FIG1_INPUT_DIR / "deepcre_random100_explanation.png").convert("RGB")
    im_mut_full = Image.open(FIG1_INPUT_DIR / "deepcre_random100_mutation.png").convert("RGB")
    im_mut_res = Image.open(FIG1_INPUT_DIR / "tutorial_canvas_10.png").convert("RGB")

    # 1. Sidebar crop from im_pred_full: x=10..485, y=355..1710 (keeps full REFERENCE ORGANISM header)
    im_sidebar = im_pred_full.crop((10, 355, 485, 1710))
    im_sidebar.save(OUT_DIR / "web_sidebar_inputs.png", dpi=(300, 300))

    # 2. Predictions crop (tabs + summary metrics + Gene Predictions Table + Genes per Chromosome):
    # x=630..3080, y=580..1950
    im_pred = im_pred_full.crop((630, 580, 3080, 1950))
    im_pred.save(OUT_DIR / "web_tab_predictions.png", dpi=(300, 300))

    # 3. Explanation crop (tabs + Average Saliency Map + Saliency Score vs Probability):
    # x=630..3080, y=580..1580
    im_expl = im_expl_full.crop((630, 580, 3080, 1580))
    im_expl.save(OUT_DIR / "web_tab_explanation.png", dpi=(300, 300))

    # 4. Mutation controls crop (tabs + Analysis Mode manual/VCF + Region + Coordinate Range + Sequence):
    im_mut_ctrl = im_mut_full.crop((630, 580, 3080, 2150))
    im_mut_ctrl.save(OUT_DIR / "web_tab_mutation_controls.png", dpi=(300, 300))
    im_mut_res.save(OUT_DIR / "web_tab_mutation_results.png", dpi=(300, 300))

    def draw_badge(x, y, num_str):
        r = 15
        draw.ellipse([x, y, x + 2 * r, y + 2 * r], fill=STEP_BG)
        tw, th = draw.textsize(num_str, font=f_step)
        draw.text((x + r - tw / 2, y + r - th / 2 - 1), num_str, fill=(255, 255, 255), font=f_step)

    # --- STEP 1: INPUT (Left column, y=135..1060) ---
    draw_badge(22, 136, "1")
    draw.text((60, 137), "Input (Sidebar)", fill=DARK, font=f_sec)

    # Place sidebar screenshot at x=24, y=180, width=215 px
    sb_w = 215
    sb_h = int(im_sidebar.height * (sb_w / im_sidebar.width))
    sb_img = add_card_border(im_sidebar.resize((sb_w, sb_h), resample=Image.BICUBIC), pad=4)
    canvas.paste(sb_img, (24, 180))

    # Callout annotations right next to the three sidebar sections
    # 1a. Reference Organism (y ~ 225)
    y1 = 215
    draw.text((256, y1), "Select Reference", fill=DARK, font=f_body)
    draw.text((256, y1 + 25), "Organism", fill=PURPLE, font=f_sub)
    draw.text((256, y1 + 54), "e.g. A. thaliana", fill=GREY, font=f_small)
    draw.text((256, y1 + 76), "(or upload FASTA/GFF)", fill=GREY, font=f_small)

    # 1b. Gene List (y ~ 405)
    y2 = 405
    draw.text((256, y2), "Upload Gene List", fill=DARK, font=f_body)
    draw.text((256, y2 + 25), "Gene IDs (.csv/.txt)", fill=PURPLE, font=f_sub)
    draw.text((256, y2 + 54), "e.g. AT4G39070,", fill=GREY, font=f_small)
    draw.text((256, y2 + 76), "AT1G72190 (or sample", fill=GREY, font=f_small)
    draw.text((256, y2 + 98), "100 random genes)", fill=GREY, font=f_small)

    # 1c. deepCRE model (y ~ 630)
    y3 = 630
    draw.text((256, y3), "Select Model", fill=DARK, font=f_body)
    draw.text((256, y3 + 25), "deepCRE model", fill=PURPLE, font=f_sub)
    draw.text((256, y3 + 54), "• Single-species", fill=GREY, font=f_small)
    draw.text((256, y3 + 76), "  (A. thaliana leaf/root)", fill=GREY, font=f_small)
    draw.text((256, y3 + 98), "• or Multi-species", fill=GREY, font=f_small)

    # Downward arrow from Step 1 to Step 2
    ax_x = 135
    draw.line([(ax_x, 825), (ax_x, 875)], fill=PURPLE, width=4)
    draw.polygon([(ax_x - 10, 870), (ax_x + 10, 870), (ax_x, 888)], fill=PURPLE)

    # --- STEP 2: PROCESSING / MODELLING (Left column, y=900..1145) ---
    draw.rectangle([22, 900, 445, 1145], fill=LIGHT_BG, outline=(205, 195, 225), width=2)
    draw_badge(36, 916, "2")
    draw.text((74, 917), "On Web Server", fill=DARK, font=f_sec)
    draw.text((38, 960), "Processing & Modelling", fill=PURPLE, font=f_sub)
    draw.text((38, 994), "• Extract flanking sequences", fill=DARK, font=f_small)
    draw.text((38, 1018), "  (gUR, gTUR, gTDR, gDR)", fill=GREY, font=f_small)
    draw.text((38, 1044), "• Predict expression class P(high)", fill=DARK, font=f_small)
    draw.text((38, 1070), "• Compute DeepSHAP saliency", fill=DARK, font=f_small)
    draw.text((38, 1096), "• Evaluate in silico mutations (ΔP)", fill=DARK, font=f_small)

    # Arrow from Step 2 to Step 3 (Outputs column on the right)
    draw.line([(445, 1020), (492, 1020)], fill=PURPLE, width=4)
    draw.polygon([(487, 1010), (487, 1030), (504, 1020)], fill=PURPLE)

    # --- STEP 3: OUTPUTS (Right column, x=510..1105, y=135..1450) ---
    out_x = 510
    out_w = 585
    draw_badge(out_x, 136, "3")
    draw.text((out_x + 38, 137), "Interactive Outputs", fill=DARK, font=f_sec)

    # 3a. Predictions tab
    y_cur = 180
    draw.text((out_x, y_cur), "Predictions", fill=PURPLE, font=f_sub)
    draw.text((out_x + 142, y_cur + 5), "— expression probability table & chromosome summary", fill=GREY, font=f_small)
    y_cur += 32
    p_h = int(im_pred.height * (out_w / im_pred.width))
    p_card = add_card_border(im_pred.resize((out_w, p_h), resample=Image.BICUBIC), pad=3)
    canvas.paste(p_card, (out_x, y_cur))
    y_cur += p_card.height + 20

    # 3b. Explanation tab
    draw.text((out_x, y_cur), "Explanation", fill=PURPLE, font=f_sub)
    draw.text((out_x + 148, y_cur + 5), "— nucleotide-resolution saliency around TSS & TTS", fill=GREY, font=f_small)
    y_cur += 32
    e_h = int(im_expl.height * (out_w / im_expl.width))
    e_card = add_card_border(im_expl.resize((out_w, e_h), resample=Image.BICUBIC), pad=3)
    canvas.paste(e_card, (out_x, y_cur))
    y_cur += e_card.height + 20

    # 3c. Mutation tab (controls + mutated vs original result)
    draw.text((out_x, y_cur), "Mutation", fill=PURPLE, font=f_sub)
    draw.text((out_x + 115, y_cur + 5), "— manual sequence edit or VCF upload & ΔP effect", fill=GREY, font=f_small)
    y_cur += 32
    tab_strip = im_mut_full.crop((630, 580, 3080, 700))
    mode_strip = im_mut_full.crop((630, 775, 3080, 1085))
    reg_strip = im_mut_full.crop((630, 1615, 3080, 2125))
    mut_comp = Image.new("RGB", (tab_strip.width, tab_strip.height + mode_strip.height + reg_strip.height), (255, 255, 255))
    mut_comp.paste(tab_strip, (0, 0))
    mut_comp.paste(mode_strip, (0, tab_strip.height))
    mut_comp.paste(reg_strip, (0, tab_strip.height + mode_strip.height))

    m_h = int(mut_comp.height * (out_w / mut_comp.width))
    m_card = add_card_border(mut_comp.resize((out_w, m_h), resample=Image.BICUBIC), pad=3)
    canvas.paste(m_card, (out_x, y_cur))
    y_cur += m_card.height + 8

    # Also place the Mutation Result (bar plot + differential saliency + delta P) right below!
    mr_h = int(im_mut_res.height * (out_w / im_mut_res.width))
    mr_card = add_card_border(im_mut_res.resize((out_w, mr_h), resample=Image.BICUBIC), pad=3)
    canvas.paste(mr_card, (out_x, y_cur))

    canvas.save(OUT_DIR / "Figure1_panel_a.png", dpi=(300, 300))
    return canvas


def main():
    im_b = make_clean_panel_b()
    im_a = make_updated_panel_a(target_height=im_b.height)

    # Combine Panel a and Panel b side-by-side with a clean subtle vertical divider
    gap = 30
    total_w = im_a.width + gap + im_b.width
    total_h = max(im_a.height, im_b.height)
    full = Image.new("RGB", (total_w, total_h), (255, 255, 255))
    full.paste(im_a, (0, 0))
    full.paste(im_b, (im_a.width + gap, 0))

    draw = ImageDraw.Draw(full)
    f_panel = ImageFont.truetype(FONT_BOLD, 50)
    f_title = ImageFont.truetype(FONT_BOLD, 30)
    # Add panel letter 'b' and title at the top of Panel b
    bx = im_a.width + gap + 10
    draw.text((bx, 12), "b", fill=(0, 0, 0), font=f_panel)
    draw.text((bx + 4, 82), "deepCRE sequence encoding & CNN architecture", fill=(35, 31, 32), font=f_title)

    # Subtle vertical separator line between Panel a and Panel b
    sep_x = im_a.width + gap // 2
    draw.line([(sep_x, 80), (sep_x, total_h - 40)], fill=(220, 220, 225), width=2)

    full.save(OUT_DIR / "Figure1.png", dpi=(300, 300))

    # Also save as 183 mm wide vector-wrapped PDF
    MM = 1.0 / 25.4
    fw_in = 183.0 * MM
    fh_in = fw_in * (total_h / total_w)
    fig = plt.figure(figsize=(fw_in, fh_in), dpi=300)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(np.array(full))
    ax.axis("off")
    fig.savefig(OUT_DIR / "Figure1.pdf", dpi=300)
    plt.close(fig)
    print(f"Saved updated Figure1.png and Figure1.pdf ({total_w}x{total_h} px)")


if __name__ == "__main__":
    main()
