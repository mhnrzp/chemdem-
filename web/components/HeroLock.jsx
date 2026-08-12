"use client";

import { useEffect, useRef } from "react";
import { heroMarkup } from "./heroMarkup";

/**
 * HeroLock renders the original navbar + hero + stats bar exactly as in
 * landing/index.html, and re-wires the original vanilla-JS behaviour
 * (mobile menu, sliding nav pill, results scroll-reveal) so it stays
 * byte-for-byte identical. Do not refactor this into JSX — the whole point
 * is pixel fidelity with the locked hero.
 */
export default function HeroLock() {
  const ref = useRef(null);

  useEffect(() => {
    const root = ref.current;
    if (!root) return;

    const byId = (id) => root.querySelector("#" + id) || document.getElementById(id);

    // --- original global handlers (inline onclick attrs resolve against window) ---
    window.toggleMobileMenu = function () {
      const menu = byId("mobileMenu");
      const btn = byId("hamburgerBtn");
      const overlay = byId("menuOverlay");
      if (!menu) return;
      if (menu.classList.contains("open")) {
        window.closeMobileMenu();
      } else {
        menu.classList.add("open");
        btn && btn.classList.add("open");
        overlay && overlay.classList.add("open");
        document.body.style.overflow = "hidden";
      }
    };
    window.closeMobileMenu = function () {
      const menu = byId("mobileMenu");
      const btn = byId("hamburgerBtn");
      const overlay = byId("menuOverlay");
      menu && menu.classList.remove("open");
      btn && btn.classList.remove("open");
      overlay && overlay.classList.remove("open");
      document.body.style.overflow = "";
    };
    const movePill = (el) => {
      const container = byId("navPill");
      const pill = byId("slidingPill");
      if (!container || !pill || !el) return;
      const containerRect = container.getBoundingClientRect();
      const elRect = el.getBoundingClientRect();
      pill.style.left = elRect.left - containerRect.left + "px";
      pill.style.width = elRect.width + "px";
    };
    window.movePill = movePill;
    window.setActive = function (el) {
      root.querySelectorAll(".nav-link").forEach((l) => l.classList.remove("active"));
      el.classList.add("active");
      movePill(el);
    };

    // --- initial pill position (original load handler) ---
    const active = root.querySelector(".nav-link.active");
    if (active) {
      const pill = byId("slidingPill");
      if (pill) {
        pill.style.transition = "none";
        movePill(active);
        requestAnimationFrame(() => {
          pill.style.transition = "";
        });
      }
    }

    // --- results scroll-reveal observer (original) ---
    const observer = new IntersectionObserver(
      (entries) => {
        entries.forEach((entry) => {
          if (entry.isIntersecting) {
            entry.target
              .querySelectorAll(".results-reveal, .success-pill-anim, .yield-fill-anim")
              .forEach((el) => el.classList.add("in-view"));
            entry.target.classList.add("in-view");
            observer.unobserve(entry.target);
          }
        });
      },
      { threshold: 0.4 }
    );
    const resultsCard = root.querySelector(".results-reveal");
    if (resultsCard) observer.observe(resultsCard);

    // --- load SmilesDrawer and render correct squaric acid chemistry ---
    const renderMols = () => {
      const mols = [
        { id: 'hero-mol-1', smiles: 'CCOC1=C(OCC)C(=O)C1=O', w: 120, h: 110 }, // DES
        { id: 'hero-mol-2', smiles: 'Nc1ccccc1',              w: 100, h: 110 }, // aniline
        { id: 'hero-mol-3', smiles: 'CCOC1=C(Nc2ccccc2)C(=O)C1=O', w: 160, h: 110 }, // product
      ];
      mols.forEach(({ id, smiles, w, h }) => {
        const canvas = root.querySelector('#' + id);
        if (!canvas) return;
        const drawer = new window.SmilesDrawer.Drawer({ width: w, height: h, bondThickness: 1.2, shortBondLength: 0.8 });
        window.SmilesDrawer.parse(smiles, tree => drawer.draw(tree, canvas, 'light'), () => {});
      });
    };

    if (window.SmilesDrawer) {
      renderMols();
    } else {
      const s = document.createElement('script');
      s.src = 'https://unpkg.com/smiles-drawer@2.0.1/dist/smiles-drawer.min.js';
      s.onload = renderMols;
      document.head.appendChild(s);
    }

    // keep pill aligned on resize
    const onResize = () => {
      const a = root.querySelector(".nav-link.active");
      if (a) movePill(a);
    };
    window.addEventListener("resize", onResize);

    return () => {
      observer.disconnect();
      window.removeEventListener("resize", onResize);
    };
  }, []);

  return <div ref={ref} dangerouslySetInnerHTML={{ __html: heroMarkup }} />;
}
