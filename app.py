# -*- coding: utf-8 -*-

import streamlit as st
import pandas as pd
import requests
from Bio import Entrez
from datetime import datetime
import time
import re
from typing import Dict, List, Optional, Tuple
import json
from bs4 import BeautifulSoup
from fake_useragent import UserAgent
import urllib.parse

# ==================== 配置与初始化 ====================
st.set_page_config(
    page_title="细胞系构建智能评估系统",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# 安全密钥配置
try:
    NCBI_EMAIL = st.secrets["NCBI_EMAIL"]
    NCBI_API_KEY = st.secrets.get("NCBI_API_KEY", "")
    APP_PASSWORD = st.secrets.get("APP_PASSWORD", "")
except Exception:
    st.error("⚠️ 请先配置 Secrets（NCBI_EMAIL 等）")
    st.stop()

# 密码保护
if APP_PASSWORD:
    if 'authenticated' not in st.session_state:
        st.session_state.authenticated = False
    if not st.session_state.authenticated:
        pwd = st.text_input("🔒 请输入访问密码", type="password")
        if pwd == APP_PASSWORD:
            st.session_state.authenticated = True
            st.rerun()
        elif pwd:
            st.error("密码错误")
        st.stop()

# 初始化 session state
if 'analysis_results' not in st.session_state:
    st.session_state.analysis_results = None
if 'search_history' not in st.session_state:
    st.session_state.search_history = []

# ==================== 核心数据获取模块 ====================

class AddgeneScraper:
    """优化的 Addgene 质粒爬取器"""
    
    def __init__(self):
        self.base_url = "https://www.addgene.org"
        self.ua = UserAgent()
        self.session = requests.Session()
        self.session.headers.update({
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
            'Accept-Language': 'en-US,en;q=0.5',
        })
    
    @st.cache_data(ttl=86400, show_spinner=False)
    def search_plasmids(_self, gene_symbol: str, max_results: int = 5) -> List[Dict]:
        """深度爬取 Addgene 质粒信息"""
        try:
            query = urllib.parse.quote(f"{gene_symbol}")
            search_url = f"{_self.base_url}/search/?q={query}&type=plasmid"
            headers = {'User-Agent': _self.ua.random}
            time.sleep(1)
            
            response = _self.session.get(search_url, headers=headers, timeout=15)
            response.raise_for_status()
            soup = BeautifulSoup(response.text, 'lxml')
            plasmids = []
            
            # 多策略解析
            result_items = (soup.find_all('article', class_='addgene-search-result') or 
                           soup.find_all('div', class_='search-result-item') or
                           soup.select('.plasmid-item'))
            
            if not result_items:
                result_items = soup.select('[data-testid="plasmid-card"]')
            
            for item in result_items[:max_results]:
                try:
                    plasmid_data = _self._parse_plasmid_card(item, gene_symbol)
                    if plasmid_data:
                        plasmids.append(plasmid_data)
                except Exception:
                    continue
            
            # 备用方案：直接访问基因页面
            if not plasmids:
                plasmids = _self._search_by_gene_page(gene_symbol)
            
            return plasmids
            
        except Exception as e:
            st.error(f"Addgene 爬取错误: {str(e)}")
            return []
    
    def _parse_plasmid_card(self, card, gene_symbol: str) -> Optional[Dict]:
        """解析单个质粒卡片"""
        try:
            link_tag = card.find('a', href=re.compile(r'/\d{5,6}/')) or card.find('a', href=True)
            if not link_tag:
                return None
            
            href = link_tag.get('href', '')
            plasmid_id_match = re.search(r'/(\d{5,6})/', href)
            if not plasmid_id_match:
                return None
            
            plasmid_id = plasmid_id_match.group(1)
            full_url = f"{self.base_url}{href}" if href.startswith('/') else href
            
            name_tag = card.find('h3') or card.find('h2') or card.find('a', class_='title')
            name = name_tag.get_text(strip=True) if name_tag else "Unknown"
            
            # 相关性检查
            if gene_symbol.lower() not in name.lower():
                desc = card.get_text().lower()
                if gene_symbol.lower() not in desc:
                    return None
            
            metadata = {
                'plasmid_id': plasmid_id,
                'name': name,
                'url': full_url,
                'insert_gene': gene_symbol,
            }
            
            # 提取表达系统
            expr_match = re.search(r'(Lenti|Retro|AAV|Mammalian|Bacterial|Insect)', 
                                  card.get_text(), re.I)
            if expr_match:
                metadata['expression_system'] = expr_match.group(1)
            
            # 提取荧光标记
            fluo_match = re.search(r'(GFP|RFP|mCherry|YFP|Luciferase|Flag|HA)', 
                                  card.get_text(), re.I)
            if fluo_match:
                metadata['tag'] = fluo_match.group(1)
            
            return metadata
            
        except Exception:
            return None
    
    def _search_by_gene_page(self, gene_symbol: str) -> List[Dict]:
        """通过基因专用页面搜索"""
        try:
            gene_url = f"{self.base_url}/browse/gene/{gene_symbol}/"
            headers = {'User-Agent': self.ua.random}
            response = self.session.get(gene_url, headers=headers, timeout=10)
            
            if response.status_code != 200:
                return []
            
            soup = BeautifulSoup(response.text, 'lxml')
            plasmids = []
            seen_ids = set()
            
            links = soup.find_all('a', href=re.compile(r'/\d{5,6}/'))
            for link in links[:5]:
                href = link.get('href', '')
                match = re.search(r'/(\d{5,6})/', href)
                if match:
                    pid = match.group(1)
                    if pid not in seen_ids:
                        seen_ids.add(pid)
                        plasmids.append({
                            'plasmid_id': pid,
                            'name': link.get_text(strip=True) or f"{gene_symbol} related",
                            'url': f"{self.base_url}/{pid}/",
                            'insert_gene': gene_symbol,
                            'source': 'Gene page'
                        })
            return plasmids
        except Exception:
            return []


class HumanAtlasScraper:
    """Human Protein Atlas 抗体信息爬取"""
    
    def __init__(self):
        self.base_url = "https://www.proteinatlas.org"
        self.ua = UserAgent()
        self.session = requests.Session()
    
    @st.cache_data(ttl=86400, show_spinner=False)
    def get_antibodies(_self, gene_symbol: str) -> List[Dict]:
        """获取经验证的抗体列表"""
        try:
            search_url = f"{_self.base_url}/{gene_symbol}"
            time.sleep(1)
            headers = {'User-Agent': _self.ua.random}
            
            response = _self.session.get(search_url, headers=headers, timeout=15)
            if response.status_code != 200:
                return []
            
            soup = BeautifulSoup(response.text, 'lxml')
            antibodies = []
            
            # 策略1：查找抗体表格
            ab_links = soup.find_all('a', href=re.compile(r'/ENSG\d+-\w+/antibody'))
            
            for link in ab_links[:5]:
                try:
                    ab_url = link.get('href')
                    if ab_url.startswith('/'):
                        ab_url = f"{_self.base_url}{ab_url}"
                    
                    ab_info = _self._parse_antibody_page(ab_url, gene_symbol)
                    if ab_info:
                        antibodies.append(ab_info)
                except Exception:
                    continue
            
            # 策略2：解析页面上的抗体列表
            if not antibodies:
                rows = soup.find_all('tr', class_=lambda x: x and 'antibody' in x.lower())
                for row in rows[:5]:
                    cells = row.find_all('td')
                    if len(cells) >= 3:
                        ab_id = cells[0].get_text(strip=True)
                        apps = _self._parse_applications(cells[1].get_text())
                        reliability = cells[2].get_text(strip=True)
                        
                        antibodies.append({
                            'gene': gene_symbol,
                            'antibody_id': ab_id,
                            'source': 'HPA',
                            'applications': apps,
                            'reliability': reliability,
                            'vendor': 'Atlas Antibodies',
                            'link': search_url
                        })
            
            return antibodies
            
        except Exception as e:
            st.warning(f"Human Atlas 查询失败: {e}")
            return []
    
    def _parse_antibody_page(self, url: str, gene_symbol: str) -> Optional[Dict]:
        """解析单个抗体页面"""
        try:
            time.sleep(0.5)
            headers = {'User-Agent': self.ua.random}
            response = self.session.get(url, headers=headers, timeout=10)
            soup = BeautifulSoup(response.text, 'lxml')
            
            # 提取抗体ID
            ab_id_match = re.search(r'(HPA\d{6}|CAB\d{6})', url)
            ab_id = ab_id_match.group(1) if ab_id_match else "Unknown"
            
            # 提取应用
            app_text = soup.get_text()
            applications = self._parse_applications(app_text)
            
            # 可靠性判断
            reliability = "Enhanced" if "enhanced" in app_text.lower() else \
                         "Supported" if "supported" in app_text.lower() else "Uncertain"
            
            return {
                'gene': gene_symbol,
                'antibody_id': ab_id,
                'source': 'HPA',
                'applications': applications,
                'reliability': reliability,
                'vendor': 'Atlas Antibodies',
                'link': url
            }
        except Exception:
            return None
    
    def _parse_applications(self, text: str) -> str:
        """解析应用类型"""
        apps = []
        text_lower = text.lower()
        app_map = {
            'western blot': 'WB', 'wb': 'WB',
            'immunohistochemistry': 'IHC', 'ihc': 'IHC',
            'immunofluorescence': 'IF', 'if': 'IF',
            'icc-if': 'ICC-IF', 'flow cytometry': 'FACS', 'facs': 'FACS'
        }
        for key, value in app_map.items():
            if key in text_lower:
                apps.append(value)
        return ', '.join(list(set(apps))) if apps else '未标注'


class BioDataFetcher:
    """生物数据获取主类"""
    
    def __init__(self, email: str, api_key: str = ""):
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        self.uniprot_base = "https://rest.uniprot.org/uniprotkb/search.json"
        self.headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'}
        self.addgene_scraper = AddgeneScraper()
        self.hpa_scraper = HumanAtlasScraper()
    
    def get_ncbi_gene_info(self, gene_symbol: str, species: str) -> Dict:
        """获取 NCBI Gene 信息"""
        try:
            term = f"{gene_symbol}[Gene Name] AND {species}[Organism]"
            handle = Entrez.esearch(db="gene", term=term, retmax=1)
            record = Entrez.read(handle)
            handle.close()
            
            if not record["IdList"]:
                return {"status": "not_found", "error": f"未找到 {gene_symbol} ({species})"}
            
            gene_id = record["IdList"][0]
            handle = Entrez.efetch(db="gene", id=gene_id, rettype="xml")
            gene_data = Entrez.read(handle)
            handle.close()
            
            gene_entry = gene_data[0]
            summary = gene_entry.get("Entrezgene_summary", "")
            
            # 致死性分析
            lethal_keywords = ["essential", "lethal", "required for cell viability", 
                             "knockout mice die", "embryonic lethal"]
            phenotype = "非必需"
            if any(kw in summary.lower() for kw in lethal_keywords):
                phenotype = "必需（潜在致死风险）"
            
            return {
                "gene_id": gene_id,
                "symbol": gene_symbol,
                "species": species,
                "description": summary[:500] if summary else "无描述",
                "phenotype": phenotype,
                "chromosome": gene_entry.get("Entrezgene_location", [{}])[0].get("Gene-location", {}).get("Gene-location_chromosome", "N/A"),
                "status": "success"
            }
        except Exception as e:
            return {"status": "error", "error": str(e)}
    
    def get_uniprot_info(self, gene_symbol: str, species: str) -> Dict:
        """获取 UniProt 蛋白信息"""
        try:
            species_map = {"Homo sapiens": "human", "Mus musculus": "mouse", "Rattus norvegicus": "rat"}
            org = species_map.get(species, species.lower())
            
            query = f"gene:{gene_symbol}+organism:{org}"
            params = {
                "query": query,
                "fields": "accession,gene_names,length,cc_subcellular_location,sequence",
                "format": "json",
                "size": 1
            }
            
            response = requests.get(self.uniprot_base, params=params, headers=self.headers, timeout=10)
            data = response.json()
            
            if not data.get("results"):
                return {"status": "not_found", "error": f"UniProt 未找到 {gene_symbol}"}
            
            protein = data["results"][0]
            accession = protein.get("primaryAccession", "")
            seq_length = protein.get("sequence", {}).get("length", 0)
            
            # 提取亚细胞定位
            loc_text = ""
            comments = protein.get("comments", [])
            for comment in comments:
                if comment.get("commentType") == "SUBCELLULAR LOCATION":
                    locations = comment.get("subcellularLocations", [])
                    locs = [loc.get("location", {}).get("value", "") for loc in locations]
                    loc_text = "; ".join([l for l in locs if l])
            
            cds_length = seq_length * 3 if seq_length else 0
            
            return {
                "uniprot_id": accession,
                "protein_length": seq_length,
                "cds_length_bp": cds_length,
                "subcellular_location": loc_text or "未标注",
                "status": "success"
            }
        except Exception as e:
            return {"status": "error", "error": str(e)}
    
    def search_pubmed_construct(self, gene_symbol: str, cell_line: Optional[str] = None, 
                              construct_type: Optional[str] = None) -> Tuple[List[Dict], int]:
        """检索 PubMed 细胞系构建文献"""
        try:
            base_query = f'{gene_symbol}[Title/Abstract]'
            
            if cell_line:
                query = f'{base_query} AND {cell_line}[Title/Abstract]'
            else:
                type_map = {
                    "overexpression": "overexpression OR over-expression OR ectopic",
                    "knockdown": "knockdown OR siRNA OR shRNA",
                    "knockout": "knockout OR CRISPR OR knock-out"
                }
                if construct_type:
                    query = f'{base_query} AND ({type_map.get(construct_type, construct_type)})'
                else:
                    query = base_query
            
            query += ' AND (cell line OR cell-line)'
            
            handle = Entrez.esearch(db="pubmed", term=query, retmax=100, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            
            pmids = record["IdList"]
            total_count = len(pmids)
            
            if not pmids:
                return [], 0
            
            fetch_ids = pmids[:10]
            handle = Entrez.efetch(db="pubmed", id=fetch_ids, rettype="abstract", retmode="xml")
            articles = Entrez.read(handle)
            handle.close()
            
            results = []
            for article in articles.get("PubmedArticle", []):
                try:
                    medline = article["MedlineCitation"]
                    article_data = medline["Article"]
                    title = article_data.get("ArticleTitle", "N/A")
                    pmid = medline.get("PMID", "N/A")
                    
                    abstract = ""
                    if "Abstract" in article_data and "AbstractText" in article_data["Abstract"]:
                        abstract = str(article_data["Abstract"]["AbstractText"])
                    
                    methods = []
                    method_keywords = ["lentiviral", "transfection", "electroporation", 
                                     "transduction", "lipofectamine", "infection"]
                    for kw in method_keywords:
                        if kw in abstract.lower():
                            methods.append(kw)
                    
                    results.append({
                        "pmid": str(pmid),
                        "title": title,
                        "methods": ", ".join(methods) if methods else "未明确提及",
                        "abstract_snippet": abstract[:200] + "..." if len(abstract) > 200 else abstract
                    })
                except Exception:
                    continue
            
            return results, total_count
        except Exception as e:
            st.error(f"PubMed 检索错误: {e}")
            return [], 0


class ConstructAnalyzer:
    """细胞系构建分析主类"""
    
    def __init__(self):
        self.fetcher = BioDataFetcher(NCBI_EMAIL, NCBI_API_KEY)
    
    def analyze_gene(self, gene_symbol: str, species: str, 
                    cell_line: Optional[str] = None, 
                    cell_species: Optional[str] = None) -> Dict:
        """执行完整分析流程"""
        
        with st.spinner(f"🔍 正在深度分析 {gene_symbol}..."):
            
            # 1. NCBI Gene
            st.text("检索 NCBI Gene...")
            ncbi_info = self.fetcher.get_ncbi_gene_info(gene_symbol, species)
            time.sleep(0.5)
            
            # 2. UniProt
            st.text("检索 UniProt...")
            uniprot_info = self.fetcher.get_uniprot_info(gene_symbol, species)
            time.sleep(0.5)
            
            # 3. Addgene（优化爬取）
            st.text("深度检索 Addgene...")
            addgene_plasmids = self.fetcher.addgene_scraper.search_plasmids(gene_symbol)
            time.sleep(0.5)
            
            # 4. HPA 抗体（优化爬取，仅人类）
            antibodies = []
            if species == "Homo sapiens":
                st.text("检索 Human Protein Atlas...")
                antibodies = self.fetcher.hpa_scraper.get_antibodies(gene_symbol)
            
            # 5. 慢病毒评估
            lentiviral = self._assess_lentiviral(ncbi_info, uniprot_info, addgene_plasmids)
            
            # 6. 文献检索
            st.text("检索细胞系构建文献...")
            literature = self._search_all_constructs(gene_symbol, cell_line)
            
            # 整合结果
            result = {
                "input_info": {
                    "gene_symbol": gene_symbol,
                    "species": species,
                    "cell_line": cell_line or "未指定",
                    "cell_species": cell_species or "未指定",
                    "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M")
                },
                "gene_function": ncbi_info,
                "protein_data": uniprot_info,
                "addgene_plasmids": addgene_plasmids,
                "antibodies": antibodies,
                "lentiviral_assessment": lentiviral,
                "cell_line_constructs": literature,
                "database_record": self._format_database_record(
                    gene_symbol, species, cell_line, ncbi_info, uniprot_info, 
                    lentiviral, literature, addgene_plasmids, antibodies
                )
            }
            
            return result
    
   def _assess_lentiviral(self, ncbi_info: Dict, uniprot_info: Dict, plasmids: List) -> Dict:
    """评估慢病毒适用性"""
    warnings = []
    recommendations = []
    score = 100  # 初始满分
    cds_len = uniprot_info.get("cds_length_bp", 0)
    
    # 致死性判断
    if ncbi_info.get("phenotype") == "必需（潜在致死风险）":
        warnings.append("⚠️ 必需基因：建议使用诱导型系统（Tet-on/off）")
        score -= 50
    
    # 序列长度判断
    if cds_len > 9000:
        warnings.append(f"⚠️ CDS {cds_len}bp 接近包装极限（10kb），包装效率可能降低")
        score -= 20
    elif cds_len > 12000:
        warnings.append(f"❌ CDS {cds_len}bp 超出慢病毒包装能力")
        score -= 80
    
    # Addgene资源
    if plasmids:
        recommendations.append(f"✅ Addgene 提供 {len(plasmids)} 个质粒")
    else:
        recommendations.append("ℹ️ Addgene 无现成质粒，需自行构建")
    
    # 根据score判断评级
    if score >= 75:
        rating = "✅ 推荐"
        suitable = True
    elif score >= 50:
        rating = "⚠️ 谨慎"
        suitable = True
    else:
        rating = "❌ 不推荐"
        suitable = False
    
    return {
        "suitable": suitable,
        "score": score,
        "warnings": warnings,
        "recommendations": recommendations,
        "overall_assessment": rating,
        "cds_length": cds_len
    }
    
    def _search_all_constructs(self, gene_symbol: str, cell_line: Optional[str]) -> Dict:
        """检索所有构建方式"""
        results = {}
        
        # 特定细胞研究
        if cell_line:
            articles, count = self.fetcher.search_pubmed_construct(gene_symbol, cell_line=cell_line)
            results["specific_cell"] = {
                "found": count > 0,
                "articles": articles,
                "total_count": count,
                "message": f"找到 {count} 篇相关文献" if count > 0 else "无完全一致的细胞系研究"
            }
        else:
            results["specific_cell"] = {"found": False, "message": "未输入细胞名称"}
        
        # 三种构建方式
        for ctype in ["overexpression", "knockdown", "knockout"]:
            articles, count = self.fetcher.search_pubmed_construct(gene_symbol, construct_type=ctype)
            results[ctype] = {
                "articles": articles[:10],
                "total_count": count,
                "methods_summary": self._extract_methods_summary(articles)
            }
        
        return results
    
    def _extract_methods_summary(self, articles: List[Dict]) -> List[str]:
        """提取常用方法"""
        all_methods = []
        for art in articles:
            methods = art.get("methods", "")
            if methods:
                all_methods.extend([m.strip() for m in methods.split(",")])
        
        from collections import Counter
        method_counts = Counter(all_methods)
        return [f"{k} ({v})" for k, v in method_counts.most_common(3)]
    
    def _format_database_record(self, gene, species, cell_line, ncbi, uniprot, 
                               lentiviral, literature, addgene_data, antibody_data) -> Dict:
        """生成数据库记录"""
        plasmid_summary = ""
        if addgene_data:
            plasmid_list = [f"{p['plasmid_id']}({p.get('expression_system', 'N/A')})" 
                           for p in addgene_data[:3]]
            plasmid_summary = "; ".join(plasmid_list)
        
        antibody_summary = ""
        if antibody_data:
            ab_list = [f"{ab['antibody_id']}({ab.get('applications', 'N/A')})" 
                      for ab in antibody_data[:3]]
            antibody_summary = "; ".join(ab_list)
        
        return {
            "基因名": gene,
            "物种": species,
            "细胞系": cell_line or "NA",
            "NCBI_Gene_ID": ncbi.get("gene_id", "NA"),
            "UniProt_ID": uniprot.get("uniprot_id", "NA"),
            "CDS长度": uniprot.get("cds_length_bp", 0),
            "基因必需性": ncbi.get("phenotype", "未知"),
            "蛋白定位": uniprot.get("subcellular_location", "未知"),
            "慢病毒适用性": lentiviral.get("overall_assessment", "未知"),
            "Addgene质粒数": len(addgene_data),
            "Addgene质粒详情": plasmid_summary,
            "可用抗体数": len(antibody_data),
            "抗体详情": antibody_summary,
            "过表达文献数": literature.get("overexpression", {}).get("total_count", 0),
            "敲低文献数": literature.get("knockdown", {}).get("total_count", 0),
            "敲除文献数": literature.get("knockout", {}).get("total_count", 0),
            "特定细胞文献": "有" if literature.get("specific_cell", {}).get("found") else "无",
            "检索日期": datetime.now().strftime("%Y-%m-%d")
        }


# ==================== Streamlit 界面 ====================

def main():
    st.title("🔬 细胞系构建智能评估系统")
    st.markdown("整合 NCBI Gene | UniProt | Addgene | Human Protein Atlas")
    
    analyzer = ConstructAnalyzer()
    
    # ==================== 第一模块：用户输入 ====================
    with st.sidebar:
        st.header("📝 第一模块：输入参数")
        
        with st.form("input_form"):
            st.subheader("基因信息（必填）")
            gene_symbol = st.text_input("基因名（如 TP53）", "").strip().upper()
            species = st.selectbox("基因物种", 
                                 ["Homo sapiens", "Mus musculus", "Rattus norvegicus"],
                                 index=0)
            
            st.subheader("细胞信息（可选）")
            cell_line = st.text_input("细胞系名称（如 HEK293T）", "").strip()
            cell_species = st.selectbox("细胞物种", 
                                       ["未指定", "Homo sapiens", "Mus musculus"], 
                                       index=0)
            
            submitted = st.form_submit_button("🔍 开始分析", use_container_width=True)
    
    # ==================== 分析执行 ====================
    if submitted:
        if not gene_symbol:
            st.error("请输入基因名")
            return
        
        result = analyzer.analyze_gene(
            gene_symbol, species, 
            cell_line if cell_line else None,
            cell_species if cell_species != "未指定" else None
        )
        
        st.session_state.analysis_results = result
        st.session_state.search_history.append(result["database_record"])
        
        st.success("✅ 分析完成！")
        st.balloons()
    
    # ==================== 展示结果 ====================
    if st.session_state.analysis_results:
        result = st.session_state.analysis_results
        
        # 基础信息
        st.divider()
        cols = st.columns(4)
        cols[0].metric("基因", result["input_info"]["gene_symbol"])
        cols[1].metric("物种", result["input_info"]["species"])
        cols[2].metric("细胞系", result["input_info"]["cell_line"])
        cols[3].metric("分析时间", result["input_info"]["analysis_date"])
        
        # ==================== 第二模块：基因与蛋白功能 ====================
        st.divider()
        st.header("🧬 第二模块：基因功能与蛋白功能")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("基因功能评估")
            gene_data = result["gene_function"]
            
            if gene_data.get("status") == "success":
                st.markdown(f"**NCBI Gene ID:** {gene_data['gene_id']}")
                st.markdown(f"**染色体:** {gene_data['chromosome']}")
                
                phenotype = gene_data.get("phenotype", "未知")
                if "必需" in phenotype:
                    st.error(f"⚠️ **致死性:** {phenotype}")
                    st.warning("建议：使用诱导型表达系统")
                else:
                    st.success(f"✅ **致死性:** {phenotype}")
                
                with st.expander("查看基因功能描述"):
                    st.info(gene_data.get("description", "无描述"))
            else:
                st.error(gene_data.get("error"))
            
            st.subheader("Addgene 质粒资源")
            plasmids = result["addgene_plasmids"]
            if plasmids:
                for p in plasmids:
                    tag_info = f" [{p.get('tag', '')}]" if p.get('tag') else ""
                    expr_info = f" ({p.get('expression_system', 'N/A')})"
                    st.markdown(f"• [{p['plasmid_id']}]({p['url']}): {p['name'][:50]}{expr_info}{tag_info}")
            else:
                st.info("未在 Addgene 找到相关质粒")
        
        with col2:
            st.subheader("蛋白信息（UniProt）")
            prot_data = result["protein_data"]
            
            if prot_data.get("status") == "success":
                st.markdown(f"**UniProt ID:** {prot_data['uniprot_id']}")
                st.markdown(f"**蛋白长度:** {prot_data['protein_length']} aa")
                st.markdown(f"**CDS长度:** {prot_data['cds_length_bp']} bp")
                st.markdown("**亚细胞定位:**")
                st.info(prot_data.get("subcellular_location", "未标注"))
            else:
                st.error(prot_data.get("error", "查询失败"))
            
            st.subheader("慢病毒包装评估")
            lentiviral = result["lentiviral_assessment"]
            
            if lentiviral["suitable"]:
                st.success(lentiviral["overall_assessment"])
            else:
                st.error(lentiviral["overall_assessment"])
            
            for w in lentiviral.get("warnings", []):
                st.warning(w)
            for r in lentiviral.get("recommendations", []):
                st.info(r)
        
        # 抗体表格
        if result["antibodies"]:
            st.subheader("经验证抗体推荐（Human Protein Atlas）")
            ab_df = pd.DataFrame(result["antibodies"])
            st.dataframe(ab_df[['antibody_id', 'applications', 'reliability', 'vendor', 'link']], 
                        use_container_width=True)
        
        # ==================== 第三模块：细胞系构建 ====================
        st.divider()
        st.header("🧫 第三模块：细胞系构建文献")
        
        lit_data = result["cell_line_constructs"]
        
        st.subheader("特定细胞系研究")
        specific = lit_data["specific_cell"]
        if specific.get("found"):
            st.success(f"✅ {specific['message']}")
            if specific.get("articles"):
                df_specific = pd.DataFrame(specific["articles"])
                st.dataframe(df_specific[['pmid', 'title', 'methods']], use_container_width=True)
        else:
            st.info(specific.get("message", "无完全一致的研究"))
        
        tabs = st.tabs(["过表达 (Overexpression)", "敲低 (Knockdown)", "敲除 (Knockout)"])
        construct_types = ["overexpression", "knockdown", "knockout"]
        
        for tab, ctype in zip(tabs, construct_types):
            with tab:
                data = lit_data[ctype]
                total = data.get("total_count", 0)
                st.markdown(f"**总文献数:** {total} 篇")
                
                if data.get("methods_summary"):
                    st.markdown(f"**常用方法:** {', '.join(data['methods_summary'])}")
                
                if data.get("articles"):
                    for article in data["articles"]:
                        with st.expander(f"{article['pmid']}: {article['title'][:60]}..."):
                            st.markdown(f"**方法:** {article['methods']}")
                            st.markdown(f"**摘要:** {article['abstract_snippet']}")
                            st.markdown(f"[PubMed链接](https://pubmed.ncbi.nlm.nih.gov/{article['pmid']}/)")
                else:
                    st.info("未找到相关文献")
        
        # ==================== 第四模块：数据库 ====================
        st.divider()
        st.header("🗄️ 第四模块：数据库整理")
        
        if st.session_state.search_history:
            db_df = pd.DataFrame(st.session_state.search_history)
            st.dataframe(db_df, use_container_width=True)
            
            col1, col2 = st.columns(2)
            with col1:
                csv = db_df.to_csv(index=False).encode('utf-8')
                st.download_button("📥 导出 CSV 数据库", csv,
                                f"cell_line_db_{datetime.now().strftime('%Y%m%d')}.csv", "text/csv")
            with col2:
                json_str = json.dumps(st.session_state.search_history, ensure_ascii=False, indent=2)
                st.download_button("📥 导出 JSON", json_str,
                                f"db_{datetime.now().strftime('%Y%m%d')}.json", "application/json")
            
            if st.button("🗑️ 清空当前会话记录"):
                st.session_state.search_history = []
                st.rerun()

if __name__ == "__main__":

    main()

